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

import edu.berkeley.cs.amplab.adam.avro.ADAMRecord
import scala.collection.mutable.ArrayBuffer

object KmerBuilder {

  /**
   * Create L-K+1 substrings of length K from a read of length L.
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
   * Create disjoint k-mers from all L-K+1 substrings of a read.
   */
  def disjointKmersFromRead (r: ADAMRecord, k: Int): ArrayBuffer[Kmer] = {
    val kmers = stringsFromRead(r, k).map(s => new Kmer(r, s))
    kmers
  }

}
