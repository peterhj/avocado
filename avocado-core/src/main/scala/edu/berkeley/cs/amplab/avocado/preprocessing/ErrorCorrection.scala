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

package edu.berkeley.cs.amplab.avocado.preprocessing

import scala.Math
//import scala.collection.JavaConversions._
import scala.collection.mutable.ArrayBuffer
import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.spark.SparkContext._
import org.apache.spark.rdd.RDD
import edu.berkeley.cs.amplab.adam.avro.ADAMRecord
import edu.berkeley.cs.amplab.avocado.assembly.Kmer

/**
 * Preprocessing stage for k-mer spectrum error correction.
 * An overview of the algorithm:
 *
 * 0. Choose a set of k-mer sizes.
 * 1. For each k, build a histogram of k-mer histogram multiplicities.
 * 2. Find the unique local minimum of each histogram and mark k-mers below it
 *    as "weak".
 * 3. Mark reads containing "weak" k-mers.
 */
object ErrorCorrection extends PreprocessingStage {

  val stageName = "errorCorrection"

  /**
   * Main method signature for a PreprocessingStage.
   *
   * @param[rdd]    an RDD of reads
   * @param[config] a SubnodeConfiguration
   * @return        an RDD of reads, with EC targets marked
   */
  def apply (rdd: RDD[ADAMRecord], config: SubnodeConfiguration): RDD[ADAMRecord] = {
    // Choose a list of k-mer sizes, and build k-mer multiplicity histograms.
    val kmer_sizes = ArrayBuffer(24)
    val kmer_mults_per_k = kmer_sizes.map(k =>
      Kmer.disjointKmersFromRDD(rdd, k)
        .map(kmer => (kmer, kmer.mult))
    )
    val mult_hists = kmer_mults_per_k.map(kmer_mults =>
      kmer_mults.map(x => (x._2, 1.asInstanceOf[Long]))
                .reduceByKey(_ + _)
                .collect
    )

    // Find the local minimum in each histogram. This is the breakpoint b/w
    // "strong" vs "weak" k-mers (to use ALLPATHS language).
    // XXX(peter, 2/14) We probably don't need to smooth the histogram,
    // especially on realistic dataset sizes. But keep smoothing in mind.
    val mult_breaks = mult_hists.map(h => {
      var mult_break: Int = -1
      val (mults, counts) = h.unzip
      for (i <- 1 until h.length - 1) {
        if (counts(i) < counts(i-1) && counts(i) < counts(i+1) && mult_break == -1) {
          mult_break = mults(i)
        }
      }
      mult_break
    })

    // Weak k-mers have multiplicity _strictly_ less than the breakpoint.
    val weak_kmers_per_k = (kmer_mults_per_k zip mult_breaks).map(pair => {
      val kmer_mults = pair._1
      val mult_break = pair._2
      val weak_kmers = kmer_mults.filter(x => x._2 < mult_break)
                                 .map(x => x._1)
      weak_kmers
    })

    // Weak reads have at least one weak k-mer, and are pointed to by their
    // k-mers. Mark these weak reads.
    weak_kmers_per_k.map(weak_kmers => {
      weak_kmers.map(kmer => {
        val k = kmer.size
        var valid_reads = (kmer.reads zip kmer.offsets).filter(x => !x._1.getFailedErrorCorrection)
        valid_reads.map(x => x._1.setRequiresErrorCorrection(true))

        /*for (r <- kmer.reads) {
          if (!r.failedErrorCorrection) {
            r.requiresErrorCorrection = true
            if (r.ecKmerOffsets == null || r.ecKmerSizes == null) {
              r.ecKmerOffsets = new ArrayBuffer[java.lang.Integer]
              r.ecKmerSizes = new ArrayBuffer[java.lang.Integer]
            }
            r.ecKmerOffsets += kmer.initialOffset
            r.ecKmerSizes += kmer.size
          }
        }*/
      })
    })

    // Editing reads is hopeless because of the sheer amount of unique k-mers.
    // Just return the reads, with some of them marked for future correction
    // such as by consensus or by realignment.
    rdd
  }

}
