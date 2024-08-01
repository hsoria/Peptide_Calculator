import streamlit as st
import pandas as pd

def calculate_amino_acid_masses(peptide_sequence, scale, resin_loading, mode, amino_acid_weights, amino_acid_comments):
    Mw_HCTU = 413.69  # g/mol
    Mw_HBTU = 379.247  # g/mol
    Mw_HOBt = 135.12
    Mw_DIEA = 129  # g/mol
    density_DIEA = 0.742  # g/mL

    if mode == "Temperature":
        amino_acid_masses = {}
        amino_acid_counts = {}
        for amino_acid in peptide_sequence:
            if amino_acid not in amino_acid_masses:
                amino_acid_masses[amino_acid] = 0
                amino_acid_counts[amino_acid] = 0
            amino_acid_masses[amino_acid] += amino_acid_weights[amino_acid] * scale * 3
            amino_acid_counts[amino_acid] += 1

        peptide_length = len(peptide_sequence)

        for amino_acid in amino_acid_masses:
            amino_acid_masses[amino_acid] /= amino_acid_counts[amino_acid]

        HCTU_pc = Mw_HCTU * 3 * scale
        amino_acid_masses["HCTU_per_coupling"] = HCTU_pc
        amino_acid_masses["resin_mg"] = (scale / resin_loading) * 1000
        amino_acid_masses["DIEA_per_coupling_µL"] = (scale * 6 * Mw_DIEA / density_DIEA)

        df = pd.DataFrame({"Quantities": amino_acid_masses})
        df = round(df, 0)
        df.rename(index=amino_acid_comments, inplace=True)

    else:
        amino_acid_masses = {}
        amino_acid_counts = {}
        for amino_acid in peptide_sequence:
            if (amino_acid not in amino_acid_masses):
                amino_acid_masses[amino_acid] = 0
                amino_acid_counts[amino_acid] = 0
            amino_acid_masses[amino_acid] += amino_acid_weights[amino_acid] * scale * 3
            amino_acid_counts[amino_acid] += 1

        peptide_length = len(peptide_sequence)

        for amino_acid in amino_acid_masses:
            amino_acid_masses[amino_acid] /= amino_acid_counts[amino_acid]

        HBTU_pc = Mw_HBTU * 3 * scale
        amino_acid_masses["HBTU_per_coupling"] = HBTU_pc

        HOBt_pc = Mw_HOBt * 3 * scale
        amino_acid_masses["HOBt_per_coupling"] = HOBt_pc

        amino_acid_masses["resin_mg"] = (scale / resin_loading) * 1000
        amino_acid_masses["DIEA_per_coupling_µL"] = (scale * 6 * Mw_DIEA / density_DIEA)
        df = pd.DataFrame({"Quantities": amino_acid_masses})
        df = round(df, 0)
        df.rename(index=amino_acid_comments, inplace=True)

    return df

def streamlit_main():
    st.title("Peptide calculator")
    st.subheader("You can choose between different high-temperature or room temperature synthesis")

    sequence = st.text_input("Enter sequence in one-letter code", value="FRGRGRG")
    scale = st.number_input("Enter scale of the synthesis [mmol]", format="%.3f", value=0.25)
    resin_loading = st.number_input("Enter scale resin loading [mmol/g]", format="%.3f", value=0.74)
    mode = st.selectbox("Synthesis method", ["RT", "Temperature"])

    add_new_amino_acid = st.selectbox("Do you want to add a new amino acid?", ["No", "Yes"])

    new_aa, new_aa_weight, new_aa_comment, new_aa_name = None, None, None, None

    if add_new_amino_acid == "Yes":
        new_aa = st.text_input("Enter new amino acid one-letter code")
        new_aa_weight = st.number_input("Enter new amino acid molecular weight", format="%.2f")
        new_aa_comment = st.text_input("Enter new amino acid comment")
        new_aa_name = st.text_input("Enter new amino acid name")

        # Validation and adding the new amino acid
        if new_aa and len(new_aa) == 1 and new_aa_weight > 0 and new_aa_comment and new_aa_name:
            st.success("New amino acid added successfully.")
        else:
            st.error("Please enter a valid one-letter code, positive molecular weight, comment, and name.")

    return sequence, scale, resin_loading, mode, new_aa, new_aa_weight, new_aa_comment, new_aa_name

def display_amino_acid_weights(amino_acid_weights, amino_acid_comments, amino_acid_names):
    amino_acids = list(amino_acid_weights.keys())
    weights = list(amino_acid_weights.values())
    comments = [amino_acid_comments[aa] for aa in amino_acid_weights]
    names = [amino_acid_names[aa] for aa in amino_acid_weights]

    df = pd.DataFrame({"Amino Acid": amino_acids, "Name": names, "Molecular weight": weights, "Comments": comments})

    centered_html = f"<div style='width: 100%; text-align: center;'>{df.to_html(index=False)}</div>"
    st.markdown(centered_html, unsafe_allow_html=True)

if __name__ == "__main__":
    amino_acid_weights = {
        "A": 311.3, "R": 648.8, "N": 354.4, "D": 411.5, "C": 585.7, "c": 391.4, 
        "E": 425.5, "Q": 368.4, "G": 297.3, "H": 619.7, "I": 353.4, "L": 353.4, 
        "K": 468.5, "M": 371.5, "F": 387.4, "P": 337.4, "S": 383.4, "T": 397.5, 
        "W": 526.58, "Y": 459.6, "y": 417.5, "V": 117.15, "Z": 101.07, "J": 178.15, 
        "U": 446.95, "O": 385.42, "o": 381.19
    }

    amino_acid_comments = {
        "A": "Fmoc-Ala", "R": "Fmoc-Arg(Pbf)", "N": "Fmoc-Asn", "D": "Fmoc-Asp(OtBu)", 
        "C": "Fmoc-Cys(Trt)", "c": "Fmoc-Cysteic", "E": "Fmoc-Glu(OtBu)", "Q": "Fmoc-Gln", 
        "G": "Fmoc-Gly", "H": "Fmoc-His(Trt)", "I": "Fmoc-Ile", "L": "Fmoc-Leu", 
        "K": "Fmoc-Lys(Boc)", "M": "Fmoc-Met", "F": "Fmoc-Phe", "P": "Fmoc-Pro", 
        "S": "Fmoc-Ser(tBu)", "T": "Fmoc-Thr(tBu)", "W": "Fmoc-Trp(Boc)", "Y": "Fmoc-Tyr(OtBu)", 
        "y": "Fmoc-Tyr(OMe)", "V": "Fmoc-Val", "Z": "Azido acetic-OH", "J": "Picolyl azide-OH", 
        "U": "Fmoc-Lys(Me3+Cl-)", "O": "Fmoc-O2Oc-OH", "o": "Fmoc-Aoa-OH"
    }

    amino_acid_names = {
        "A": "Alanine", "R": "Arginine", "N": "Asparagine", "D": "Aspartic acid", 
        "C": "Cysteine", "c": "Cysteic acid", "E": "Glutamic acid", "Q": "Glutamine", 
        "G": "Glycine", "H": "Histidine", "I": "Isoleucine", "L": "Leucine", 
        "K": "Lysine", "M": "Methionine", "F": "Phenylalanine", "P": "Proline", 
        "S": "Serine", "T": "Threonine", "W": "Tryptophan", "Y": "Tyrosine", 
        "y": "Methyl tyrosine", "V": "Valine", "Z": "Azido acetic acid", "J": "Picolyl azide", 
        "U": "Trimethyllysine", "O": "Ornithine", "o": "Aminocyclobutanecarboxylic acid"
    }

    sequence, scale, resin_loading, mode, new_aa, new_aa_weight, new_aa_comment, new_aa_name = streamlit_main()

    if new_aa and len(new_aa) == 1 and new_aa_weight > 0 and new_aa_comment and new_aa_name:
        amino_acid_weights[new_aa] = new_aa_weight
        amino_acid_comments[new_aa] = new_aa_comment
        amino_acid_names[new_aa] = new_aa_name

    df = calculate_amino_acid_masses(sequence, scale, resin_loading, mode, amino_acid_weights, amino_acid_comments)
    st.write(df)

    st.title("Amino Acid Legend")
    display_amino_acid_weights(amino_acid_weights, amino_acid_comments, amino_acid_names)
