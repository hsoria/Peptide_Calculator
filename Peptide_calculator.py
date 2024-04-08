import streamlit as st
import pandas as pd
import numpy as np

import pandas as pd




def calculate_amino_acid_masses(peptide_sequence, scale, resin_loading, mode):
    amino_acid_weights = {
        "A": 311.3,
        "R": 648.8, #Fmoc-Arg(pbf)
        "N": 354.4, #Fmoc-Asn
        "D": 411.5, #Fmoc-Asp(OtBu)
        "C": 585.7,#Fmoc-Cys(Trt)
        "c": 391.4, #Fmoc-Cysteic
        "E": 425.5, #Fmoc-Glu(OtBu)
        "Q": 368.4, #Fmoc-Gln
        "G": 297.3,
        "H": 619.7, #Fmoc-HisTrt
        "I": 353.4, 
        "L": 353.4,
        "K": 468.5, #Fmoc-Lys(Boc)
        "M": 371.5, #FMoc-Met
        "F": 387.4,
        "P": 337.4,
        "S": 383.4, #Fmoc-Ser(tBu)
        "T": 397.5, #Fmoc-Thr(tBu)
        "W": 426.5, #Fmoc-Trp 
        "Y": 459.6,
        "y": 417.5,
        "V": 117.15,
        "Z": 101.07, #Azidoacetic acid
        "u": 446.95,
        "O": 385.42
    }
    Mw_HCTU = 413.69 #g/mol
    Mw_HBTU = 379.247 #g/mol
    Mw_HOBt = 135.12
    Mw_DIEA = 129 #g/mol
    density_DIEA = 0.742 #g/mL

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
        HCTU_mg = (peptide_length+1)* HCTU_pc
        amino_acid_masses["HCTU_per_coupling"] =  HCTU_pc
        amino_acid_masses["resin_mg"] =  (scale/resin_loading)*1000
        amino_acid_masses["DIEA_per_coupling_µL"] = (scale * 6 * Mw_DIEA / density_DIEA)

        df = pd.DataFrame({"Quantities": amino_acid_masses})
        df = round(df, 0)

        
    else: 
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
            
        HBTU_pc = Mw_HBTU * 4 * scale
        HBTU_mg = (peptide_length+1)* HBTU_pc
        amino_acid_masses["HBTU_per_coupling"] =  HBTU_pc

        HOBt_pc = Mw_HOBt * 4 * scale
        HOBt_mg = (peptide_length+1)* HOBt_pc
        amino_acid_masses["HOBt_per_coupling"] =  HOBt_pc

        amino_acid_masses["resin_mg"] =  (scale/resin_loading)*1000
        amino_acid_masses["DIEA_per_coupling_µL"] = (scale * 6 * Mw_DIEA / density_DIEA)
        df = pd.DataFrame({"Quantities": amino_acid_masses})
        df = round(df, 0)

    
    return df


def streamlit_main():
    st.title("Peptide calculator")
    st.subheader("You can choose between different high-temperature or room temperature synthesis")

    sequence = st.text_input("Enter sequence in one-letter code", value = "FRGRGRG")
    scale = st.number_input("Enter scale of the synthesis [mmol]", format="%.3f", value = 0.25)
    resin_loading = st.number_input("Enter scale resin loading [mmol/g]", format="%.3f", value = 0.74)
    mode = st.selectbox("Synthesis method", ["RT", "Temperature"])

    return sequence, scale, resin_loading, mode


if __name__ == "__main__":
    sequence, scale, resin_loading, mode = streamlit_main()
    df = calculate_amino_acid_masses(sequence, scale, resin_loading, mode)
    
    st.write(df)


def display_amino_acid_weights(amino_acid_weights):
    amino_acids = list(amino_acid_weights.keys())
    weights = list(amino_acid_weights.values())
    comments = ["Fmoc-Ala", "Fmoc-Arg(Pbf)", "Fmoc-Asn", "Fmoc-Asp(OtBu)", "Fmoc-Cys(Trt)", "Fmoc-Cysteic", "Fmoc-Glu(OtBu)",
                "Fmoc-Gln", "Fmoc-Gly", "Fmoc-His(Trt)", "Fmoc-Ile", "Fmoc-Leu", "Fmoc-Lys(Boc)", "Fmoc-Met", "Fmoc-Phe", "Fmoc-Pro",
                "Fmoc-Ser(tBu)", "Fmoc-Thr(tBu)", "Fmoc-Trp(Boc)", "Fmoc-Tyr(OtBu)", "Fmoc-Tyr(OMe)",
                  "Fmoc-Val", "Azido acetic-OH", "Picolyl azide-OH", "Fmoc-Lys(Me3+Cl-)", "Fmoc-O2Oc-OH"]
    
    df = pd.DataFrame({"Amino Acid": amino_acids, "Molecular weight": weights, "Molecule": comments})
    
    # Adjust the table to the window width and center it
    centered_html = f"<div style='width: 100%; text-align: center;'>{df.to_html(index=False)}</div>"
    st.markdown(centered_html, unsafe_allow_html=True)

if __name__ == "__main__":
    amino_acid_weights = {
        "A": 311.3, #
        "R": 648.8, #
        "N": 354.4, #Fmoc-Asn
        "D": 411.5, #Fmoc-Asp(OtBu)
        "C": 585.7,#Fmoc-Cys(Trt)
        "c": 391.4, #Fmoc-Cysteic
        "E": 425.5, #Fmoc-Glu(OtBu)
        "Q": 368.4, #Fmoc-Gln
        "G": 297.3, #Fmoc-Gly
        "H": 619.7, #Fmoc-HisTrt
        "I": 353.4, #Fmoc-Ile
        "L": 353.4, #Fmoc-Leu
        "K": 468.5, #Fmoc-Lys(Boc)
        "M": 371.5, #Fmoc-Met
        "F": 387.4, #Fmoc-Phe
        "P": 337.4, #Fmoc-Pro
        "S": 383.4, #Fmoc-Ser(tBu)
        "T": 397.5, #Fmoc-Thr(tBu)
        "W": 526.58, #Fmoc-Trp 
        "Y": 459.6, #Fmoc-Tyr()
        "y": 417.5, #Fmoc-Tyr(OMe)
        "V": 117.15, #Fmoc-Val
        "Z": 101.07, #Azidoacetic acid
        "z":178.15, #Picolyl azide
        "u":446.95, 
        "O":385.42, 

    }
    st.title("Amino Acid legend")
    display_amino_acid_weights(amino_acid_weights)
