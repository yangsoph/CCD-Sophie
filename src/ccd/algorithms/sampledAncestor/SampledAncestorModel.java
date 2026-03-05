package ccd.algorithms.sampledAncestor;

import ccd.algorithms.regularisation.CCDRegularisationStrategy;

public enum SampledAncestorModel {

    ExtendedClade("ExtendedClade"),
    ExtendedLeaf("ExtendedLeaf"),
    NotApplicable("NotApplicable");

    private final String name;

    SampledAncestorModel(String name) {
        this.name = name;
    }

    @Override
    public String toString() {
        return this.name;
    }
}