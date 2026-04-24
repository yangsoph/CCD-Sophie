package ccd.model;

import ccd.model.bitsets.BitSet;

public class CladeSA extends Clade {

    public CladeSA(BitSet cladeInBits, AbstractCCD abstractCCD) {
        super(cladeInBits, abstractCCD);
    }

    public CladeSA(AbstractCCD abstractCCD) {
        super(abstractCCD);
    }

}
