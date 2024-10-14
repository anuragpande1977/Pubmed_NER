import subprocess
import pandas as pd
import streamlit as st
from Bio import Entrez, Medline
from io import BytesIO
import matplotlib.pyplot as plt
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from transformers import AutoTokenizer, AutoModelForTokenClassification, pipeline
from collections import Counter

# Load the Hugging Face model for biomedical NER
def load_huggingface_model():
    tokenizer = AutoTokenizer.from_pretrained("d4data/biomedical-ner-all")
    model = AutoModelForTokenClassification.from_pretrained("d4data/biomedical-ner-all")
    return pipeline("ner", model=model, tokenizer=tokenizer)

ner_model = load_huggingface_model()

# Extract entities using the Hugging Face model
def extract_entities(text):
    entities = ner_model(text)
    return [ent['word'] for ent in entities if 'DISEASE' in ent['entity']]

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

# Create an Excel file
def save_to_excel(articles):
    output = BytesIO()
    data = [{'Title': art.get('TI', 'N/A'), 'Abstract': art.get('AB', 'N/A')} for art in articles]
    pd.DataFrame(data).to_excel(output, index=False)
    output.seek(0)
    return output

# Streamlit UI
st.title("PubMed Biomedical NER Search")
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
