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
import edu.berkeley.cs.amplab.adam.avro.ADAMRecord

class Kmer (read: ADAMRecord, val string: String) {

  var reads = new ArrayBuffer[ADAMRecord](1)

  reads += read

  override def equals(that: Any): Boolean = {
    that match {
      case kmer: Kmer => string == kmer.string
      case _ => false
    }
  }

  def mult(): Int = {
    reads.length
  }

}
