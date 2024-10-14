import subprocess
import pandas as pd
import streamlit as st
from Bio import Entrez, Medline
from io import BytesIO
import matplotlib.pyplot as plt
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import spacy
from collections import Counter

# Function to ensure the model is installed and loaded
def load_model():
    try:
        nlp = spacy.load("en_core_web_sm")
    except OSError:
        st.warning("Model not found. Installing en_core_web_sm...")
        subprocess.run(["python", "-m", "spacy", "download", "en_core_web_sm"], check=True)
        nlp = spacy.load("en_core_web_sm")
    return nlp

# Load the spaCy model
nlp = load_model()

# Function to extract entities from text using spaCy
def extract_entities(text):
    doc = nlp(text)
    return [ent.text for ent in doc.ents if ent.label_ == "DISEASE"]

# Fetch articles from PubMed
def fetch_abstracts(query, num_articles, email):
    Entrez.email = email
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

# Create an Excel file in memory
def save_to_excel(articles):
    output = BytesIO()
    data = [{'Title': art.get('TI', 'N/A'), 'Abstract': art.get('AB', 'N/A')} for art in articles]
    pd.DataFrame(data).to_excel(output, index=False)
    output.seek(0)
    return output

# Streamlit UI
st.title("PubMed NER Search")
email = st.text_input("Enter your email:")
search_term = st.text_input("Search term:")
num_articles = st.number_input("Number of articles:", 1, 100, 10)

if st.button("Fetch and Analyze"):
    if email and search_term:
        articles = fetch_abstracts(search_term, num_articles, email)
        if articles:
            st.download_button("Download Results", save_to_excel(articles), "results.xlsx")
            for article in articles:
                abstract = article.get('AB', '')
                diseases = extract_entities(abstract)
                st.write(f"Diseases: {', '.join(diseases) if diseases else 'None'}")
        else:
            st.write("No results found.")
