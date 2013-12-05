/*
 * Copyright (c) 2013. Regents of the University of California
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package edu.berkeley.cs.amplab.avocado.calls.reads

import org.apache.spark.{SparkContext, Logging}
import org.apache.spark.rdd.RDD
import edu.berkeley.cs.amplab.adam.avro.{ADAMRecord, ADAMVariant, ADAMGenotype}
import scala.collection.mutable.{ArrayBuffer, HashMap, HashSet, PriorityQueue}
import scala.math._

class HMMAlignment (_ref: CharSequence, _test: CharSequence, _ref_start: Int) {
  val ref_haplotype: CharSequence = _ref
  val test_haplotype: CharSequence = _test

  var ref_start = _ref_start

  val padded_ref_len = ref_haplotype.length + 1
  val padded_test_len = test_haplotype.length + 1

  val stride = padded_ref_len
  val mat_size = padded_test_len * padded_ref_len

  // TODO(peter, 12/4) HMM parameters from the snp_prob and indel_prob?
  //val epsilon = 0.0
  //val delta = 0.0
  //val tau = 0.0
  val eta = -log10(1.0 + test_haplotype.length) // Note: want to use the _test_
  // haplotype length here, not the ref length.
  val d = -4.0
  val e = -4.0

  // TODO(peter, 12/3) using a default quality for now, should use read scores.
  val q_phred = 40
  val q = pow(10.0, q_phred / -10.0)
  val p = 1.0 - q
  val prior_match = p
  val prior_indel = q
  val match_to_match = p
  val indel_to_match = p
  val match_to_indel = q
  val indel_to_indel = q

  var trace: ArrayBuffer[Char] = null
  var matches: ArrayBuffer[Double] = null
  var inserts: ArrayBuffer[Double] = null
  var deletes: ArrayBuffer[Double] = null

  def computeLogLikelihood(): Double = {
    // TODO(peter, 12/4) shortcut b/w global and local alignment: use a custom
    // start position in the reference haplotype.
    trace = new ArrayBuffer[Char]
    matches = new ArrayBuffer[Double]
    inserts = new ArrayBuffer[Double]
    deletes = new ArrayBuffer[Double]
    matches(ref_start) = 2.0 * eta
    inserts(ref_start) = Double.NegativeInfinity
    deletes(ref_start) = Double.NegativeInfinity
    for (i <- 0 to padded_test_len) {
      for (j <- ref_start to padded_ref_len) {
        if (i > 0 || j > 0) {
          val idx = i * stride + j

          val m_prior = 0.0 // TODO(peter, 12/4) different match prior if snp or not.
          val m_match = matches((i-1) * stride + (j-1))
          val m_insert = inserts((i-1) * stride + (j-1))
          val m_delete = deletes((i-1) * stride + (j-1))
          val m = m_prior + max(m_match, max(m_insert, m_delete))

          val ins_match = matches((i-1) * stride + j) + d
          val ins_insert = inserts((i-1) * stride + j) + e
          val ins = max(ins_match, ins_insert)

          val del_match = matches(i * stride + (j-1)) + d
          val del_delete = deletes(i * stride + (j-1)) + e
          val del = max(del_match, del_delete)

          val best = max(m, max(ins, del))
          val t = if (best == m) {
            'M'
          } else if (best == ins) {
            'I'
          } else if (best == del) {
            'D'
          } else {
            '.'
          }

          trace(idx) = t
          matches(idx) = m
          inserts(idx) = ins
          deletes(idx) = del
        }
      }
    }
    val v = max(matches(mat_size - 1), max(inserts(mat_size - 1), deletes(mat_size - 1)))
    v
  }

  def computeAlignment(): Unit = {
    // TODO(peter, 12/4) just make a cigar string?
    for (i <- 1 to padded_test_len) {
      for (j <- (ref_start + 1) to padded_ref_len) {
        val idx = i * stride + j
        val t = trace(idx)
        t match {
          case 'M' => {}
          case 'I' => {}
          case 'D' => {}
          case _ => {}
        }
      }
    }
  }
}

// For our purposes, a read is a list of kmers.
class AssemblyRead (_record: ADAMRecord) {
  val record: ADAMRecord = _record
  var kmers: ArrayBuffer[Kmer] = null
}

// A kmer prefix is a string of length k-1.
class KmerPrefix(_string: CharSequence) {
  val string: CharSequence = _string
}

// A kmer has a prefix of length k-1 and a unit length suffix.
class Kmer (_prefix: KmerPrefix, _suffix: CharSequence, _left: KmerVertex, _right: KmerVertex) {
  val prefix: KmerPrefix = _prefix
  val suffix: CharSequence = _suffix
  var left: KmerVertex = _left
  var right: KmerVertex = _right
  var reads: HashSet[AssemblyRead] = new HashSet[AssemblyRead]
  var mult: Int = 1
}

class KmerVertex {
  var left: HashSet[Kmer] = new HashSet[Kmer]
  var right: HashSet[Kmer] = new HashSet[Kmer]
}

class KmerPath (_edges: ArrayBuffer[Kmer]) {
  val edges = _edges
  var mult_sum: Int = 0
  for (e <- edges) {
    mult_sum += e.mult
  }
}

object KmerPathOrdering extends Ordering[KmerPath] {
  def compare(path1: KmerPath, path2: KmerPath): Int = {
    // Kmer paths are ordered by increasing mult sum.
    if (path1.mult_sum < path2.mult_sum) {
      return 1
    }
    else if (path1.mult_sum > path2.mult_sum) {
      return -1
    }
    0
  }
}

// A mate pair is a pair of reads.
//class AssemblyMatePair (_r1: AssemblyRead, _r2: AssemblyRead) {
//  var read1: AssemblyRead = _r1
//  var read2: AssemblyRead = _r2
//}

class KmerGraph (_klen: Int, _region_len: Int) {
  val klen: Int = _klen
  val spur_threshold: Int = klen // TODO(peter, 11/26) how to choose thresh?

  //var reads: HashMap[CharSequence,AssemblyRead] = null
  var reads: List[AssemblyRead] = null

  // Convenient to explicitly have the graph source and sink.
  val source = new KmerVertex
  val sink = new KmerVertex

  // The actual kmer graph consists of unique K-1 prefixes and kmers connected
  // by vertices.
  var prefixes = new HashMap[CharSequence,KmerPrefix]
  var kmers = new HashMap[KmerPrefix,ArrayBuffer[Kmer]]

  //var all_paths = new HashSet[KmerPath]
  var all_paths = new PriorityQueue[KmerPath]()(KmerPathOrdering)

  def insertReadKmers(r: AssemblyRead): Unit = {
    // Construct L-K+1 kmers, initially connected to the source and sink.
    val readlen = r.record.sequence.length
    val offsets = ArrayBuffer.range(0, readlen - klen + 1)
    val ks = offsets.map(idx => {
      val prefix_str = r.record.sequence.subSequence(idx, idx + klen - 1)
      val prefix = prefixes.getOrElseUpdate(prefix_str, new KmerPrefix(prefix_str))
      val suffix = r.record.sequence.subSequence(idx + klen - 1, idx + klen)
      var k = new Kmer(prefix, suffix, source, sink)
      k.reads.add(r)
      kmers.getOrElseUpdate(k.prefix, new ArrayBuffer[Kmer]) += k
      k
    })

    // Add a vertex in between each adjacent pair of kmers.
    (ks, ks.drop(1)).zipped.map((k1, k2) => {
      var v = new KmerVertex
      // TODO(peter, 11/26) check the order of k1 and k2!
      k1.right = v
      k2.left = v
      v.left.add(k1)
      v.right.add(k2)
    })
    ks.head.left = source
    ks.last.right = sink
    source.right.add(ks.head)
    sink.left.add(ks.last)

    r.kmers = ks
  }

  def insertReads(read_group: List[ADAMRecord]): Unit = {
    reads = read_group.map(x => new AssemblyRead(x))
    reads.map(r => insertReadKmers(r))
  }

  def mergeVertices(v1: KmerVertex, v2: KmerVertex): Unit = {
    // Merge v2 into v1.
    v2.left.map(k => v1.left.add(k))
    v2.right.map(k => v1.right.add(k))
  }

  def exciseVertexKmer(v: KmerVertex, k: Kmer): Unit = {
    v.left.remove(k)
    v.right.remove(k)
  }

  def exciseKmer(k: Kmer): Unit = {
  }

  def exciseVertex(v: KmerVertex): Unit = {
  }

  def connectGraph(): Unit = {
    for ((prefix, ks) <- kmers) {
      // For each prefix, each suffix has an arbitrary "canonical" kmer.
      var canon_ks = new HashMap[CharSequence,Kmer]
      for (k <- ks) {
        var canon_k = canon_ks.getOrElseUpdate(k.suffix, k)
        if (k != canon_k) {
          // Consolidate equivalent kmers together (i.e., same prefix and suffix)
          // with additive mults. Also fixup the vertices.
          mergeVertices(canon_k.left, k.left)
          mergeVertices(canon_k.right, k.right)
          exciseVertexKmer(canon_k.left, k)
          exciseVertexKmer(canon_k.right, k)
          canon_k.reads ++= k.reads
          canon_k.mult += k.mult
        }
      }
    }
  }

  def removeSpurs(): Unit = {
    // Remove all kmer paths with length <= spur_threshold, connected to the
    // graph source and sink.
    // TODO(peter, 11/27) make sure the path lengths are initialized and
    // updated correctly!
    for (sk <- source.right) {
      var path_len = 0
      var k = sk
      while (k.right.right.size == 1) {
        path_len += 1 //k.suffix.length
        k = k.right.right.head
      }
      if (path_len <= spur_threshold) {
        k = k.left.left.head
        while (k != source) {
          val lk = k.left.left.head
          exciseVertex(k.left)
          exciseKmer(k)
          k = lk
        }
      }
    }
    for (sk <- sink.left) {
      var path_len = 0
      var k = sk
      while (k.left.left.size == 1) {
        path_len += 1 //k.left.left.head.suffix.length
        k = k.left.left.head
      }
      if (path_len <= spur_threshold) {
        k = k.right.right.head
        while (k != source) {
          val rk = k.right.right.head
          exciseVertex(k.right)
          exciseKmer(k)
          k = rk
        }
      }
    }
  }

  def threadReads(): Unit = {
    // TODO(peter, 11/27)
    // Unwind "frayed rope" patterns using individual reads.
    // At the moment, done directly on the kmer graph, but could also condense
    // the kmer graph into the repeat graph, and thread on that.
    var convergings = new HashMap[AssemblyRead,ArrayBuffer[Kmer]]
    for ((_, ks) <- kmers) {
      for (k <- ks) {
        if (k.left.left.size > 1 &&
            k.left.right.size == 1)
        {
          for (r <- k.reads) {
            convergings.getOrElseUpdate(r, new ArrayBuffer[Kmer]) += k
          }
        }
      }
    }
    for ((r, ks) <- convergings) {
      for (k <- ks) {
        // TODO(peter, 11/27)
        // Find the diverging kmer, then duplicate all kmers in between, given
        // sufficient evidence.
      }
    }
  }

  //def threadMatePairs(): Unit = {
  //}

  def enumerateAllPaths(): Unit = {
    // Do DFS to enumerate and score all paths through the final graph.
    val max_depth = 400 - klen + 1
    var edges = new ArrayBuffer[Kmer]
    def allPathsDFS(v: KmerVertex, depth: Int): Unit = {
      if (v == sink) {
        val path = new KmerPath(edges)
        all_paths.enqueue(path)
      }
      else if (depth <= max_depth) {
        for (k <- v.right) {
          edges += k
          allPathsDFS(k.right, depth + 1)
          edges.remove(edges.length - 1)
        }
      }
    }
    allPathsDFS(source, 0)
  }
}

/**
 * Phase (diploid) haplotypes with kmer assembly on active regions.
 */
class ReadCallAssemblyPhaser extends ReadCall {
  override val callName = "AssemblyPhaser"

  val kmer_len = 20
  val region_len = 200

  def assemble(read_group: List[ADAMRecord]): KmerGraph = {
    var kmer_graph = new KmerGraph(kmer_len, region_len)
    kmer_graph.insertReads(read_group)
    kmer_graph.connectGraph()
    //kmer_graph.removeSpurs()
    kmer_graph.enumerateAllPaths()
    kmer_graph
  }

  def phaseAssembly(kmer_graph: KmerGraph): List[(ADAMVariant, List[ADAMGenotype])] = {
    null
  }

  /**
   * Method signature for variant calling operation.
   *
   * @param[in] pileupGroups An RDD containing reads.
   * @return An RDD containing called variants.
   */
  override def call (pileupGroups: RDD [ADAMRecord]): RDD [(ADAMVariant, List[ADAMGenotype])] = {
    null
  }
}
