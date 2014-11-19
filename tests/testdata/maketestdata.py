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