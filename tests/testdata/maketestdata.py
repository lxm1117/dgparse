# Copyright (c) 2014, Rebecca R. Murphy
# All rights reserved.

# * Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions are met:

# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.

# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.


import struct

# making damamged files to test error handling


def make_duplicates():
    segments = [0]

    for i in segments:
        fl = "pDONR223 empty vector.dna"
        with open(fl, "rb") as f:
            with open("duplicate.dna", "w") as outfile:
                segment = f.read(5)
                while segment:
                    seg, seg_len = struct.unpack('>BI', segment)
                    data = f.read(seg_len)
                    if seg != i:
                        outfile.write(segment)
                        outfile.write(data)
                    else:
                        for s in range(2):
                            outfile.write(segment)
                            outfile.write(data)
                    print seg, seg_len
                    segment = f.read(5)


def make_truncations():
    segments = [0]

    for i in segments:
        fl = "pDONR223 empty vector.dna"
        with open(fl, "rb") as f:
            with open("truncated.dna", "w") as outfile:
                segment = f.read(5)
                while segment:
                    seg, seg_len = struct.unpack('>BI', segment)
                    data = f.read(seg_len)
                    if seg != i:
                        outfile.write(segment)
                        outfile.write(data)
                    else:
                        outfile.write(segment)
                        outfile.write(data[:100])
                    print seg, seg_len
                    segment = f.read(5)


def delete_segments():
    segments = [0, 10, 5, 8, 6, 2, 9]

    for i in segments:
        fl = "pDONR223 empty vector.dna"
        with open(fl, "rb") as f:
            with open("test_no%s.dna" % i, "w") as outfile:
                segment = f.read(5)
                while segment:
                    seg, seg_len = struct.unpack('>BI', segment)
                    data = f.read(seg_len)
                    if seg != i:
                        outfile.write(segment)
                        outfile.write(data)
                    print seg, seg_len
                    segment = f.read(5)


def main():
    delete_segments()
    make_truncations()
    make_duplicates()


if __name__ == "__main__":
    main()