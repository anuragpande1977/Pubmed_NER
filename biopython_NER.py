import subprocess
import pandas as pd
import streamlit as st
from Bio import Entrez, Medline
from io import BytesIO
import matplotlib.pyplot as plt
from collections import Counter
from transformers import AutoTokenizer, AutoModelForTokenClassification, pipeline
import re

# Email validation function
def is_valid_email(email):
    return re.match(r"[^@]+@[^@]+\.[^@]+", email)

# Load a more specialized Hugging Face NER model for scientific text
@st.cache_resource
def load_huggingface_model():
    tokenizer = AutoTokenizer.from_pretrained("allenai/scibert_scivocab_cased")
    model = AutoModelForTokenClassification.from_pretrained("d4data/biomedical-ner-all")
    return pipeline("ner", model=model, tokenizer=tokenizer, aggregation_strategy="simple")

ner_model = load_huggingface_model()

# Extract disease-related entities
def extract_entities(text):
    entities = ner_model(text)
    extracted_terms = [
        entity['word'] for entity in entities 
        if any(label in entity['entity'] for label in ['DISEASE', 'B-DISEASE', 'I-DISEASE'])
    ]
    return extracted_terms

# PubMed query construction
def construct_query(search_term, mesh_term, article_type):
    query = f"({search_term}) AND {article_type}[pt]"
    if mesh_term:
        query += f" AND {mesh_term}[MeSH Terms]"
    return query

# Fetch articles from PubMed with pagination support
def fetch_abstracts(query, num_articles, email):
    Entrez.email = email
    articles = []
    try:
        batch_size = 100
        for start in range(0, num_articles, batch_size):
            handle = Entrez.esearch(
                db="pubmed", term=query, retstart=start, retmax=batch_size
            )
            result = Entrez.read(handle)
            ids = result['IdList']
            handle.close()

            if ids:
                handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
                records = Medline.parse(handle)
                articles.extend(records)
                handle.close()

        return articles

    except Exception as e:
        st.write(f"An error occurred: {e}")
        return []

# Save articles to Excel
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

# Plot top 15 disease terms by frequency with improved aesthetics
def plot_top_disease_frequency(disease_list):
    if disease_list:
        disease_freq = Counter(disease_list)
        top_15 = disease_freq.most_common(15)
        df = pd.DataFrame(top_15, columns=["Disease", "Frequency"])

        plt.figure(figsize=(10, 6))
        bars = plt.bar(
            df["Disease"], df["Frequency"], 
            color=plt.cm.viridis(df["Frequency"] / max(df["Frequency"]))
        )
        plt.xticks(rotation=45, ha='right')
        plt.xlabel("Disease Terms")
        plt.ylabel("Frequency")
        plt.title("Top 15 Disease Terms in Abstracts")
        plt.tight_layout()
        st.pyplot(plt)
    else:
        st.write("No disease terms found.")

# Streamlit User Interface
st.title("PubMed NER Search and Top 15 Disease Term Frequency")

email = st.text_input("Enter your email for PubMed access:")
search_term = st.text_input("Enter the search term:")
mesh_term = st.text_input("Enter an optional MeSH term (leave blank if not needed):")
article_type = st.selectbox("Select article type:", ["Clinical Trial", "Review", "Meta-Analysis", "Case Report"])
num_articles = st.number_input("Number of articles to fetch:", min_value=1, max_value=500, value=10)

if st.button("Search"):
    if not is_valid_email(email):
        st.write("Please enter a valid email address.")
    elif not search_term:
        st.write("Search term cannot be empty.")
    else:
        query = construct_query(search_term, mesh_term, article_type)
        with st.spinner("Fetching articles..."):
            articles = fetch_abstracts(query, num_articles, email)

        if articles:
            excel_data = save_to_excel(articles)
            st.download_button(
                label="Download Results as Excel",
                data=excel_data,
                file_name="pubmed_results.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )

            all_entities = []
            for article in articles:
                abstract = article.get('AB', '')
                entities = extract_entities(abstract)
                all_entities.extend(entities)

            if all_entities:
                st.write("### Top 15 Disease Term Frequency:")
                plot_top_disease_frequency(all_entities)
            else:
                st.write("No disease terms found.")
        else:
            st.write("No articles found.")

