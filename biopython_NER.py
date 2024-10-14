import subprocess
import pandas as pd
import streamlit as st
from Bio import Entrez, Medline
from io import BytesIO
import spacy
from collections import Counter

# Install and load the spaCy model (en_core_web_sm)
def load_model():
    try:
        nlp = spacy.load("en_core_web_sm")
    except OSError:
        st.warning("Model not found. Installing en_core_web_sm...")
        subprocess.run(["python", "-m", "spacy", "download", "en_core_web_sm"])
        nlp = spacy.load("en_core_web_sm")
    return nlp

# Load the model
nlp = load_model()

# Define article types for PubMed search
article_types = {
    "Clinical Trials": "Clinical Trial[pt]",
    "Meta-Analysis": "Meta-Analysis[pt]",
    "Randomized Controlled Trials": "Randomized Controlled Trial[pt]",
    "Reviews": "Review[pt]",
    "Systematic Reviews": "Systematic Review[pt]",
    "Case Reports": "Case Reports[pt]",
    "Observational Studies": "Observational Study[pt]",
}

# Function to construct the query with optional MeSH terms
def construct_query(search_term, mesh_term, article_type):
    query = f"({search_term}) AND {article_types[article_type]}"
    if mesh_term:
        query += f" AND {mesh_term}[MeSH Terms]"
    return query

# Fetch articles from PubMed using Entrez API
def fetch_abstracts(query, num_articles, email):
    Entrez.email = email
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=num_articles)
        result = Entrez.read(handle)
        ids = result['IdList']
        handle.close()

        if not ids:
            st.write("No articles found.")
            return []

        handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
        records = Medline.parse(handle)
        articles = list(records)
        handle.close()
        return articles

    except Exception as e:
        st.write(f"An error occurred: {e}")
        return []

# Function to extract disease-related terms using spaCy NER
def extract_disease_terms(text):
    doc = nlp(text)
    return [ent.text for ent in doc.ents if ent.label_ == "DISEASE"]

# Create an Excel file with the results
def save_to_excel(articles):
    output = BytesIO()
    data = [
        {
            "Title": article.get("TI", "No Title"),
            "Abstract": article.get("AB", "No Abstract"),
            "Authors": ", ".join(article.get("AU", "No Authors")),
            "Journal": article.get("TA", "No Journal"),
            "Publication Date": article.get("DP", "No Date")
        }
        for article in articles
    ]
    df = pd.DataFrame(data)
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        df.to_excel(writer, index=False)
    output.seek(0)
    return output

# Streamlit User Interface
st.title("PubMed Research and Disease Term Finder")

email = st.text_input("Enter your email for PubMed access:")
search_term = st.text_input("Enter the search term:")
mesh_term = st.text_input("Enter an optional MeSH term (leave blank if not needed):")
article_type = st.selectbox("Select article type:", list(article_types.keys()))
num_articles = st.number_input("Number of articles to fetch:", min_value=1, max_value=100, value=10)

if st.button("Search"):
    if email and search_term:
        query = construct_query(search_term, mesh_term, article_type)
        articles = fetch_abstracts(query, num_articles, email)

        if articles:
            # Save the articles to Excel
            excel_data = save_to_excel(articles)
            st.download_button(
                label="Download Results as Excel",
                data=excel_data,
                file_name="pubmed_results.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )

            # Display disease terms from abstracts
            st.write("### Extracted Disease Terms from Abstracts:")
            for article in articles:
                abstract = article.get('AB', '')
                disease_terms = extract_disease_terms(abstract)
                if disease_terms:
                    st.write(f"Title: {article.get('TI', 'No Title')}")
                    st.write(f"Disease Terms: {', '.join(disease_terms)}")
                else:
                    st.write(f"Title: {article.get('TI', 'No Title')}")
                    st.write("No disease terms found.")
        else:
            st.write("No articles found.")
    else:
        st.write("Please provide both email and search term.")

