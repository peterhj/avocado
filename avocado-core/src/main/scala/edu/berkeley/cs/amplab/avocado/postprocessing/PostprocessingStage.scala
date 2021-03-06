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

package edu.berkeley.cs.amplab.avocado.postprocessing

import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.spark.rdd.RDD
import edu.berkeley.cs.amplab.adam.models.ADAMVariantContext
import edu.berkeley.cs.amplab.avocado.stats.AvocadoConfigAndStats

private[postprocessing] trait PostprocessingStage {

  val stageName: String

  def apply (rdd: RDD[ADAMVariantContext], 
             stats: AvocadoConfigAndStats,
             config: SubnodeConfiguration): RDD[ADAMVariantContext]

}
