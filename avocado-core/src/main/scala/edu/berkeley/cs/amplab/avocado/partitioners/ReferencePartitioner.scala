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

package edu.berkeley.cs.amplab.avocado.partitioners

import org.apache.commons.configuration.{HierarchicalConfiguration, SubnodeConfiguration}
import org.apache.spark.rdd.RDD
import edu.berkeley.cs.amplab.adam.avro.ADAMRecord
import edu.berkeley.cs.amplab.avocado.stats.AvocadoConfigAndStats

trait ReferencePartitionerCompanion {

  val partitionerName: String

  protected def apply (rdd: RDD[ADAMRecord], 
                       subnodeConfiguration: SubnodeConfiguration,
                       stats: AvocadoConfigAndStats): ReferencePartitioner 

  final def apply (rdd: RDD[ADAMRecord], 
                   globalConfig: HierarchicalConfiguration,
                   partitionName: String,
                   stats: AvocadoConfigAndStats): ReferencePartitioner = {
    
    val config: SubnodeConfiguration = globalConfig.configurationAt(partitionName)

    apply(rdd, config, stats)
  }

}

trait ReferencePartitioner {
  
  val companion: ReferencePartitionerCompanion

  def computePartitions (): PartitionSet
  
}
