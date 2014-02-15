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

import scala.collection.mutable.ArrayBuffer
import org.apache.spark.SparkContext._
import org.apache.spark.rdd.RDD
import edu.berkeley.cs.amplab.adam.avro.ADAMRecord

class Kmer (val string: String, read: ADAMRecord, offset: Int) {

  var reads = new ArrayBuffer[ADAMRecord](1)
  var offsets = new ArrayBuffer[Int](1)

  reads += read
  offsets += offset

  override def equals (that: Any): Boolean = {
    that match {
      case kmer: Kmer => string == kmer.string
      case _ => false
    }
  }

  def size (): Int = {
    string.length
  }

  def mult (): Int = {
    reads.length
  }

  def initialRead (): ADAMRecord = {
    reads(0)
  }

  def initialOffset (): Int = {
    offsets(0)
  }

}

object Kmer {

  /**
   * Create L-K+1 substrings of length K from a read of length L.
   *
   * @param[r]  a read
   * @param[k]  the k-mer size
   * @return    an array of nonunique k-mer strings
   */
  def stringsFromRead (r: ADAMRecord, k: Int): ArrayBuffer[String] = {
    val readlen = r.sequence.length
    var strings = new ArrayBuffer[String](readlen - k + 1)
    for (i <- 0 until readlen - k + 1) {
      strings(i) = r.sequence.toString.substring(i, i + k)
    }
    strings
  }

  /**
   * Create (disconnected, unconsolidated) k-mers from all L-K+1 substrings of
   * a read.
   *
   * @param[r]  a read
   * @param[k]  the k-mer size
   * @return    an array of nonunique and disjoint k-mers
   */
  def disjointKmersFromRead (r: ADAMRecord, k: Int): ArrayBuffer[Kmer] = {
    val readlen = r.sequence.length
    var kmers = new ArrayBuffer[Kmer](readlen - k + 1)
    for (i <- 0 until readlen - k + 1) {
      val string = r.sequence.toString.substring(i, i + k)
      kmers(i) = new Kmer(string, r, i)
    }
    kmers
  }

  /**
   * Consolidate a collection of k-mers sharing a common string value.
   *
   * @param[kmers]  a collection of k-mers with the same string
   * @return        a canonical, unique k-mer replacing the previous k-mers
   */
  def consolidateEqualKmers (kmers: Seq[Kmer]): Kmer = {
    var canon_kmer = kmers(0)
    for (kmer <- kmers) {
      if (!kmer.eq(canon_kmer)) {
        canon_kmer.reads ++= kmer.reads
      }
    }
    canon_kmer
  }

  /**
   * Create (disconnected) consolidated k-mers from a collection of reads.
   *
   * @param[rdd]  a collection of reads
   * @param[k]    the k-mer size
   * @return      an array of consolidated but disjoint kmers
   */
  def disjointKmersFromReads (rs: Seq[ADAMRecord], k: Int): ArrayBuffer[Kmer] = {
    ArrayBuffer[Kmer](rs.flatMap(r => disjointKmersFromRead(r, k))
      .groupBy(kmer => kmer.string)
      .toSeq
      .map(x => consolidateEqualKmers(x._2)) : _*)
  }

  /**
   * Create (disconnected) consolidated k-mers from an RDD of reads.
   *
   * @param[rdd]  an RDD of reads
   * @param[k]    the k-mer size
   * @return      an RDD of consolidated but disjoint kmers
   */
  def disjointKmersFromRDD (rdd: RDD[ADAMRecord], k: Int): RDD[Kmer] = {
    rdd.flatMap(r => disjointKmersFromRead(r, k))
       .map(kmer => (kmer.string, kmer))
       .groupByKey
       .map(x => consolidateEqualKmers(x._2))
  }
}
