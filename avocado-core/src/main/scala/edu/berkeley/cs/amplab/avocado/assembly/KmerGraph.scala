/*
 * Copyright (c) 2014. Regents of the University of California
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

package edu.berkeley.cs.amplab.avocado.assembly

import scala.collection.mutable.{ArrayBuffer, HashMap, HashSet}
import edu.berkeley.cs.amplab.adam.avro.ADAMRecord

class KmerEdge (var kmer: Kmer) {

  var v1 = KmerVertex.leftOfEdge(this)
  var v2 = KmerVertex.rightOfEdge(this)

}

object KmerVertex {

  def leftOfEdge (edge: KmerEdge): KmerVertex = {
    var v = new KmerVertex
    v.outEdges += edge
    v
  }

  def rightOfEdge (edge: KmerEdge): KmerVertex = {
    var v = new KmerVertex
    v.inEdges += edge
    v
  }

}

class KmerVertex {

  var inEdges = new ArrayBuffer[KmerEdge]
  var outEdges = new ArrayBuffer[KmerEdge]

}

class KmerGraph (val uncorrReads: Seq[ADAMRecord], val readLength: Int, val kmerSize: Int) {

  val reads = uncorrReads.flatMap(r => SpectralErrorCorrector.correctRead(r))
  var kmers = Kmer.disjointKmersFromReads(reads, kmerSize)

  var edges = kmers.map(kmer => new KmerEdge(kmer))
  var vertexSet = new HashSet[KmerVertex]
  var kmerMap = new HashMap[String, KmerEdge]

  var readKmersMap = new HashMap[ADAMRecord, ArrayBuffer[KmerEdge]]

  init

  def init () = {
    edges.map(edge => kmerMap.put(edge.kmer.string, edge))
    threadKmers
    computeReadKmersMap
  }

  def mergeVertices (vertex: KmerVertex, other: KmerVertex): KmerVertex = {
    vertex.inEdges ++= other.inEdges
    vertex.outEdges ++= other.outEdges
    vertex
  }

  def joinEdges (e1: KmerEdge, e2: KmerEdge) = {
    val v = mergeVertices(e1.v2, e2.v1)
    e1.v2 = v
    e2.v1 = v
  }

  def threadKmers () = {
    // The only connections b/w k-mers we allow are the ones originally present
    // in the reads themselves.
    for (r <- reads) {
      for (i <- 0 until readLength - kmerSize) {
        val s1 = r.getSequence.toString.substring(i, i + kmerSize)
        val s2 = r.getSequence.toString.substring(i + 1, i + kmerSize + 1)
        val e1 = kmerMap.get(s1).get
        val e2 = kmerMap.get(s2).get
        joinEdges(e1, e2)
      }
    }
    for (edge <- edges) {
      vertexSet.add(edge.v1)
      vertexSet.add(edge.v2)
    }
  }

  def findSources (): ArrayBuffer[KmerVertex] = {
    var sources = new ArrayBuffer[KmerVertex]
    for (vertex <- vertexSet) {
      if (vertex.inEdges.length == 0) {
        sources += vertex
      }
    }
    sources
  }

  def computeReadKmersMap () = {
    for (e <- edges) {
      for (r <- e.kmer.reads) {
        readKmersMap.getOrElse(r, new ArrayBuffer[KmerEdge]) += e
      }
    }
  }

}
