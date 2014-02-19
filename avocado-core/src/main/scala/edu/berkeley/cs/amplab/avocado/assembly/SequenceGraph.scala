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

import scala.collection.mutable.{ArrayBuffer, HashMap, HashSet, Queue}
import edu.berkeley.cs.amplab.adam.avro.ADAMRecord

class Kpath {

  var kmers = new ArrayBuffer[KmerEdge]
  var edge: Option[SequenceEdge] = None

  private var string: Option[String] = None

  def length (): Int = {
    if (kmers.length == 0) {
      0
    }
    else {
      kmers(0).kmer.string.length + kmers.length - 1
    }
  }

}

class SequenceEdge (var kpath: Kpath, var v1: SequenceVertex, var v2: SequenceVertex) {

  kpath.edge = Some(this)
  v1.outEdges += this
  v2.inEdges += this

}

class SequenceVertex {

  var inEdges = new ArrayBuffer[SequenceEdge]
  var outEdges = new ArrayBuffer[SequenceEdge]

}

class SequencePath {

  var edges = new ArrayBuffer[SequenceEdge]

}

class SequenceReadPair (val r1: ADAMRecord, val r2: ADAMRecord) {

  var spanningPaths = new ArrayBuffer[SequencePath]

}

class SequenceFragment (val readPair: SequenceReadPair) {
}

class SequenceGraph (val kmerGraph: KmerGraph, val insertMean: Double, val insertSdev: Double) {

  val reads = kmerGraph.reads
  val readLength = kmerGraph.readLength
  val kmerSize = kmerGraph.kmerSize

  var kpaths = new ArrayBuffer[Kpath]
  var edges = new ArrayBuffer[SequenceEdge]
  var kmerEdgeMap = new HashMap[KmerEdge, (Kpath, Int)]
  var kmerVertexMap = new HashMap[KmerVertex, SequenceVertex]

  var unpairedReads = new ArrayBuffer[ADAMRecord]
  var readPairs = new ArrayBuffer[SequenceReadPair]
  var completeFragments = new ArrayBuffer[SequenceFragment]

  init

  def init () = {
    var visited_kmer_tips = new HashSet[KmerEdge]
    val sources = kmerGraph.findSources
    sources.map(source => makeRecursive(visited_kmer_tips, source))
    pairReads
    fillAllPaths
    fillFragments
  }

  def kpathWithInitialKmer(initialKmer: KmerEdge): (Kpath, KmerVertex) = {
    var kmer = initialKmer
    var kpath = new Kpath
    var offset = 0
    while (kmer.v2.inEdges.length == 1 && kmer.v2.outEdges.length == 1) {
      kpath.kmers += kmer
      kmerEdgeMap.put(kmer, (kpath, offset))
      kmer = kmer.v2.outEdges(0)
    }
    kpath.kmers += kmer
    kmerEdgeMap.put(kmer, (kpath, offset))
    (kpath, kmer.v2)
  }

  def makeRecursive (visited_kmer_tips: HashSet[KmerEdge], initial_kmer_vertex: KmerVertex): Unit = {
    var initial_vertex = kmerVertexMap.getOrElse(initial_kmer_vertex, new SequenceVertex)
    for (kmer <- initial_kmer_vertex.outEdges) {
      if (!visited_kmer_tips.contains(kmer)) {
        visited_kmer_tips.add(kmer)
        val (kpath, terminal_kmer_vertex) = kpathWithInitialKmer(kmer)
        var terminal_vertex = kmerVertexMap.getOrElse(terminal_kmer_vertex, new SequenceVertex)
        edges += new SequenceEdge(kpath, initial_vertex, terminal_vertex)
        makeRecursive(visited_kmer_tips, terminal_kmer_vertex)
      }
    }
  }

  def pairReads () = {
    val readsByName = reads.groupBy(r => r.getReadName.toString)
    unpairedReads ++= readsByName.filter(kv => kv._2.length == 1)
                                 .map(kv => kv._2(0))
    readPairs ++= readsByName.filter(kv => kv._2.length == 2)
                             .map(kv => new SequenceReadPair(kv._2(0), kv._2(1)))
  }

  def fillAllPaths () = {
    for (pair <- readPairs) {
      val initial_string = pair.r1.getSequence.toString.substring(0, kmerSize)
      val initial_kmer = kmerGraph.kmerMap.get(initial_string).get
      val initial_edge = kmerEdgeMap.get(initial_kmer).get._1.edge.get
      val terminal_string = pair.r2.getSequence.toString.substring(readLength - kmerSize - 1, readLength - 1)
      val terminal_kmer = kmerGraph.kmerMap.get(terminal_string).get
      val terminal_edge = kmerEdgeMap.get(terminal_kmer).get._1.edge.get

      // Find all paths between initial_vertex and terminal_vertex by doing
      // distance-bounded DFS.
      fillPathRecursive(pair, initial_edge, terminal_edge, initial_edge.kpath.length - kmerEdgeMap.get(initial_kmer).get._2)
    }
  }

  def fillPathRecursive (pair: SequenceReadPair, edge: SequenceEdge, terminal_edge: SequenceEdge, dist: Int): Unit = {
    for (e2 <- edge.v2.outEdges) {
      if (e2 == terminal_edge) {
      }
      else if (dist > Math.ceil(insertMean + 3.0 * insertSdev).toInt) {
      }
    }
  }

  def fillFragments () = {
  }

}
