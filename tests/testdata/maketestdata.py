import struct

def main():
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


if __name__ == "__main__":
    main()