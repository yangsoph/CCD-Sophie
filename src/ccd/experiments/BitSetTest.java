package ccd.experiments;

import ccd.model.bitsets.BitSet;

import java.io.IOException;

public class BitSetTest {
    public static void main(String[] args) throws IOException {
        BitSet cladeInBits = BitSet.newBitSet(10);
        // cladeInBits.set(0);
        cladeInBits.set(5);
        System.out.println("size = " + cladeInBits.size());
        System.out.println("length = " + cladeInBits.length());
        System.out.println("get(0) = " + cladeInBits.get(0));
        System.out.println("get(1) = " + cladeInBits.get(1));
        System.out.println("get(2) = " + cladeInBits.get(2));
        BitSet sub = cladeInBits.getSubset(4,6);
        System.out.println("sub length = " + sub.length());
        System.out.println(sub.get(1));
    }
}
