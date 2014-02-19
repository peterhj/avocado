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

package edu.berkeley.cs.amplab.avocado.calls.reads

import scala.collection.mutable.ArrayBuffer
import edu.berkeley.cs.amplab.adam.avro.ADAMRecord
import edu.berkeley.cs.amplab.adam.models.ADAMVariantContext
import edu.berkeley.cs.amplab.avocado.assembly.{KmerGraph, SequenceGraph}
import edu.berkeley.cs.amplab.avocado.calls.VariantCallCompanion
import edu.berkeley.cs.amplab.avocado.stats.{AvocadoConfigAndStats, InsertDist}
import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.spark.{SparkContext, Logging}
import org.apache.spark.rdd.RDD

object ReadCallAssembler extends VariantCallCompanion {

  val callName = "ReadAssembler"

  def apply (stats: AvocadoConfigAndStats,
             config: SubnodeConfiguration): ReadCallAssembler = {
    new ReadCallAssembler(stats)
  }

}

class ReadCallAssembler (val stats: AvocadoConfigAndStats) extends ReadCall {

  val companion = ReadCallAssembler

  override def isCallable () = true

  def call (pileupGroups: RDD[ADAMRecord]): RDD[ADAMVariantContext] = {
    null // FIXME
  }

  def callRegion (region: Seq[ADAMRecord]): ArrayBuffer[ADAMVariantContext] = {
    val kmer_graph = new KmerGraph(region, stats.readLength, 24)
    val seq_graph = new SequenceGraph(kmer_graph, stats.insertDist.mean, stats.insertDist.sdev)

    null // FIXME
  }

}
