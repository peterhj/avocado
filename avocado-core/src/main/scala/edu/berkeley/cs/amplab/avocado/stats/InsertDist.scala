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

package edu.berkeley.cs.amplab.avocado.stats

import edu.berkeley.cs.amplab.adam.avro.ADAMRecord
import org.apache.spark.SparkContext._
import org.apache.spark.rdd.RDD

private[stats] object InsertDist {

  def apply(rdd: RDD[ADAMRecord]): InsertDist = {
    // FIXME(peter, 2/16) proper way to test for "good" mate pair?
    val inserts = rdd.filter(r => r.getReadMapped && r.getMateMapped && r.getFirstOfPair)
                     .map(r => {
                       var ins = r.getMateAlignmentStart - r.getStart
                       ins = if (ins >= 0) { ins } else { -ins }
                       ins
                     })
    // XXX(peter) watch out for long overflow if the insert is too large.
    val inserts_count = inserts.count
    val inserts_sum = inserts.reduce(_ + _)
    val insert_mean = inserts_sum.toDouble / inserts_count.toDouble
    val inserts_square_sum = inserts.map(ins => (ins - insert_mean) * (ins - insert_mean))
                                    .reduce(_ + _)
    val insert_sdev = Math.sqrt(inserts_square_sum.toDouble / inserts_count.toDouble)
    new InsertDist(insert_mean, insert_sdev)
  }

}

class InsertDist (val mean: Double, val sdev: Double) {
  // TODO(peter) empirical dist?
}
