package starbeast2.app.beauti;

import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.StandardData;
import beast.base.evolution.datatype.UserDataType;
import beast.base.parser.PartitionContext;
import beast.base.spec.evolution.alignment.FilteredAlignment;
import beastfx.app.inputeditor.BeautiDoc;
import morphmodels.app.beauti.BeautiMorphModelAlignmentProvider;
import starbeast2.SpeciesTree;
import starbeast2.StarBeastTaxonSet;

import java.util.*;

@Description("Class for creating new partitions for morphological data to be edited by AlignmentListInputEditor")
public class StarBeastMorphModelAlignmentProvider extends BeautiMorphModelAlignmentProvider {
    @Override
    public void processAlignment(Alignment alignment, List<BEASTInterface> filteredAlignments, boolean ascertained, BeautiDoc doc) throws Exception {
        StarBeastTaxonSet ts = (StarBeastTaxonSet) doc.pluginmap.get("taxonsuperset");
        ts.alignmentInput.set(alignment);
        ts.initAndValidate();

        SpeciesTree st = (SpeciesTree) doc.pluginmap.get("Tree.t:Species");
        st.m_taxonset.set(ts);
        st.makeCaterpillar(0, 1, false);

        // super.processAlignment(alignment, filteredAlignments, ascertained, doc);
        // TODO: revert to original super.processAlignment, after MM bug is fixed.
        // https://github.com/CompEvol/morph-models/issues/21
        splitByStateSpace(alignment, filteredAlignments, ascertained, doc);
    }

    /**
     * Copy of {@link BeautiMorphModelAlignmentProvider#processAlignment} that creates
     * beast.base.spec.evolution.alignment.FilteredAlignment instead of the legacy
     * beast.base.evolution.alignment.FilteredAlignment, to match the
     * StarBEASTMKTrait subtemplate. Otherwise BEAUti's XMLParser rejects the
     * pre-registered $(n) object (error 105) and the partition is never created.
     * Remove once morph-models creates the spec type itself.
     */
    private void splitByStateSpace(Alignment alignment, List<BEASTInterface> filteredAlignments, boolean ascertained, BeautiDoc doc) throws Exception {
        Map<Integer, List<Integer>> stateSpaceMap = new HashMap<>();

        int initialSiteCount = alignment.getSiteCount();
        int maxNrOfStates = 0;

        // distinguish between StandardData and others
        if (alignment.getDataType() instanceof StandardData) {
            // determine state space size by interrogating StandardData data-type
            StandardData dataType = (StandardData) alignment.getDataType();
            for (int i = 0; i < alignment.getSiteCount(); i++) {
                int nrOfStates;
                int nrOfStatesPresented = calcNumberOfStates(alignment, i);
                if (dataType.charStateLabelsInput.get().size() > i && dataType.charStateLabelsInput.get().get(i).getStateCount() > 0) {
                    // this assumes there is a charStateLabel with the state description for this site
                    nrOfStates = dataType.charStateLabelsInput.get().get(i).getStateCount();
                    if (nrOfStatesPresented > nrOfStates) {
                        throw new Exception("The number of states in character " + (i+1) + " is larger than in " +
                                "the description. It should be less or equal.");
                    }
                } else {
                    // deal with the case where there is no charStateLabel or there is no state description
                    if (nrOfStatesPresented < 2) {
                        throw new RuntimeException("Cannot determine the number of possible states for character " +
                                (i+1) + ". \n There is no character description and there are fewer than two states for " +
                                "this character in the matrix. \n Please specify the number of possible states for " +
                                "characters in CHARSTATELABELS block");
                    }
                    nrOfStates = nrOfStatesPresented;
                }
                if (!stateSpaceMap.containsKey(nrOfStates)) {
                    stateSpaceMap.put(nrOfStates, new ArrayList<>());
                    maxNrOfStates = Math.max(maxNrOfStates, nrOfStates);
                }
                stateSpaceMap.get(nrOfStates).add(i);
            }
        } else {
            // determine state space size by counting states in each site
            for (int i = 0; i < alignment.getSiteCount(); i++) {
                int nrOfStates = calcNumberOfStates(alignment, i);
                if (!stateSpaceMap.containsKey(nrOfStates)) {
                    stateSpaceMap.put(nrOfStates, new ArrayList<>());
                    maxNrOfStates = Math.max(maxNrOfStates, nrOfStates);
                }
                stateSpaceMap.get(nrOfStates).add(i);
            }
        }

        if (ascertained) {
            StringBuilder seqToAccountForAscertainment = new StringBuilder();
            for (int i = 0; i < maxNrOfStates; i++) {
                seqToAccountForAscertainment.append(i);
            }
            for (Sequence seq : alignment.sequenceInput.get()) {
                String newSequenceValue = seq.dataInput.get() + seqToAccountForAscertainment;
                seq.dataInput.setValue(newSequenceValue, seq);
            }
            alignment.initAndValidate();
        }

        String tree = alignment.getID();
        String clock = alignment.getID();

        // create filtered alignments, one per state space size
        for (Integer nrOfStates : stateSpaceMap.keySet()) {
            String name = alignment.getID() + nrOfStates;

            // create filter range, collapsing consecutive sites into ranges
            StringBuilder range = new StringBuilder();
            List<Integer> sites = stateSpaceMap.get(nrOfStates);
            for (int i = 0; i < sites.size(); i++) {
                int site = sites.get(i);
                if (i == sites.size()-1 || sites.get(i+1) != site+1) {
                    range.append((site + 1) + ",");
                } else {
                    if (range.length() == 0 || range.charAt(range.length()-1) != '-') {
                        range.append((site + 1) + "-");
                    }
                }
            }
            if (ascertained) {
                range.append((initialSiteCount+1) + "-" + (initialSiteCount+nrOfStates) + ",");
            }
            range.deleteCharAt(range.length() - 1);

            // create data type
            DataType.Base dataType;
            if (alignment.getDataType() instanceof StandardData) {
                StandardData base = (StandardData) alignment.getDataType();
                dataType = new StandardData();
                ((StandardData) dataType).initByName("nrOfStates", nrOfStates,
                        "ambiguities", base.listOfAmbiguitiesInput.get());
            } else {
                // TODO deal with ambiguous codes
                StringBuilder codeMap = new StringBuilder();
                for (int i = 0; i < nrOfStates; i++) {
                    codeMap.append(i + " = " + i + ", ");
                }
                codeMap.append("? =");
                for (int i = 0; i < nrOfStates; i++) {
                    codeMap.append(" " + i);
                }
                UserDataType userDataType = new UserDataType();
                userDataType.initByName("states", nrOfStates, "codelength", 1, "codeMap", codeMap.toString());
                dataType = userDataType;
            }
            dataType.setID("morphDataType." + name);
            doc.addPlugin(dataType);

            FilteredAlignment data = new FilteredAlignment();
            if (ascertained) {
                data.isAscertainedInput.setValue(true, data);
                data.excludefromInput.setValue(stateSpaceMap.get(nrOfStates).size(), data);
                data.excludetoInput.setValue(stateSpaceMap.get(nrOfStates).size()+nrOfStates, data);
            }
            data.initByName("data", alignment, "filter", range.toString(), "userDataType", dataType);
            data.setID(name);
            doc.addPlugin(data);

            // link trees and clock models, and create treelikelihood for each state space
            PartitionContext context = new PartitionContext(name, name, clock, tree);
            try {
                doc.addAlignmentWithSubnet(context, template.get());
            } catch (Exception e) {
                e.printStackTrace();
            }

            filteredAlignments.add(data);
        }
    }

    /**
     * calculate number of states for a site by determining the set of unique
     * characters at that site.
     */
    private int calcNumberOfStates(Alignment alignment, int site) {
        int[] pattern = alignment.getPattern(alignment.getPatternIndex(site));
        Set<Integer> states = new HashSet<>();
        DataType dataType = alignment.getDataType();
        for (int k : pattern) {
            if (k >= 0 && !dataType.getCode(k).equals("?") && !dataType.getCode(k).equals("-")) {
                for (int m : dataType.getStatesForCode(k)) {
                    states.add(m);
                }
            }
        }
        return states.size();
    }
}
