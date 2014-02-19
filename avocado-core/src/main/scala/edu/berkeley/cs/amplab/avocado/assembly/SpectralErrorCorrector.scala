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

//import scala.collection.JavaConversions._
//import scala.collection.mutable.{ArrayBuffer, Buffer}
import edu.berkeley.cs.amplab.adam.avro.ADAMRecord

object SpectralErrorCorrector {

  def correctRead(r: ADAMRecord): List[ADAMRecord] = {
    List(r)

    /*if (!r.requiresErrorCorrection ||
        r.ecKmerOffsets == null ||
        r.ecKmerSizes == null)
    {
      return r
    }

    var ec = ADAMRecord.newBuilder(r).build
    val offsets: Buffer[java.lang.Integer] = r.ecKmerOffsets
    val sizes: Buffer[java.lang.Integer] = r.ecKmerSizes
    val pairs = (offsets zip sizes).map(x => (x._1.intValue, x._2.intValue))
    for ((offset, size) <- pairs) {
    }

    ec*/
  }

}
