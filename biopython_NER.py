import os
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

# Function to install the SciSpacy model if it's not available
def load_model():
    try:
        nlp = spacy.load("en_ner_bc5cdr_md")
    except OSError:
        st.warning("SciSpacy model not found. Installing it now...")
        subprocess.run([
            "pip", "install",
            "https://github.com/allenai/scispacy-models/releases/download/en_core_sci_sm-0.5.0/en_core_sci_sm-0.5.0.tar.gz"
        ])
        nlp = spacy.load("en_core_sci_sm")  # Use an alternative model
    return nlp

# Load the SciSpacy model
nlp = load_model()

# Define PubMed article types and search tags
article_types = {
    "Clinical Trials": "Clinical Trial[pt]",
    "Meta-Analysis": "Meta-Analysis[pt]",
    "Randomized Controlled Trials": "Randomized Controlled Trial[pt]",
    "Reviews": "Review[pt]",
    "Systematic Reviews": "Systematic Review[pt]",
    "Case Reports": "Case Reports[pt]",
    "Observational Studies": "Observational Study[pt]",
    "Comparative Studies": "Comparative Study[pt]",
    "Editorials": "Editorial[pt]",
}

# Construct query with optional MeSH term
def construct_query(search_term, mesh_term, choice):
    chosen_article_type = article_types[choice]
    query = f"({search_term}) AND {chosen_article_type}"
    if mesh_term:
        query += f" AND {mesh_term}[MeSH Terms]"
    return query

# Fetch articles from PubMed
def fetch_abstracts(query, num_articles, email):
    Entrez.email = email
    st.write(f"Searching PubMed with query: {query}")
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=num_articles)
        result = Entrez.read(handle)
        handle.close()

        ids = result['IdList']
        if not ids:
            st.write("No articles found matching the criteria.")
            return []

        handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
        records = Medline.parse(handle)
        articles = list(records)
        handle.close()
        return articles
    except Exception as e:
        st.write(f"An error occurred: {e}")
        return []

# Generate Excel file in memory
def save_to_excel(articles):
    output = BytesIO()
    data = [{
        'Title': article.get('TI', 'No title available'),
        'Authors': ', '.join(article.get('AU', 'No authors available')),
        'Abstract': article.get('AB', 'No abstract available'),
        'Publication Date': article.get('DP', 'No publication date available'),
        'Journal': article.get('TA', 'No journal available')
    } for article in articles]
    df = pd.DataFrame(data)
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        df.to_excel(writer, index=False)
    output.seek(0)
    return output

# Count articles by year
def count_articles_by_year(articles):
    year_count = {}
    for article in articles:
        pub_date = article.get('DP', 'No publication date available')
        if pub_date != 'No publication date available':
            year = pub_date.split()[0]
            if year.isdigit():
                year_count[year] = year_count.get(year, 0) + 1
    return year_count

# Plot pie chart with Plotly
def plot_publication_years_pie_chart(year_count):
    if year_count:
        fig = go.Figure(data=[go.Pie(
            labels=list(year_count.keys()),
            values=list(year_count.values()),
            hoverinfo='label+value',
            textinfo='label+value'
        )])
        fig.update_layout(title="Distribution of Articles by Publication Year", title_x=0.5)
        st.plotly_chart(fig)

# Extract disease terms from abstracts using SciSpacy
def extract_disease_terms(abstract_text):
    doc = nlp(abstract_text)
    return [ent.text for ent in doc.ents if ent.label_ == "DISEASE"]

# Analyze and plot disease term frequency
def analyze_disease_frequency(articles, top_n=10):
    disease_list = []
    for article in articles:
        abstract = article.get('AB', 'No abstract available')
        if abstract != 'No abstract available':
            disease_list.extend(extract_disease_terms(abstract))
    disease_freq = Counter(disease_list)
    if disease_freq:
        df = pd.DataFrame(disease_freq.items(), columns=["Disease", "Frequency"])
        df = df.sort_values(by="Frequency", ascending=False).head(top_n)
        plt.figure(figsize=(10, 6))
        plt.bar(df["Disease"], df["Frequency"])
        plt.xlabel('Diseases')
        plt.ylabel('Frequency')
        plt.title(f'Top {top_n} Disease Term Frequencies')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        st.pyplot(plt)
    else:
        st.write("No disease terms found in the abstracts.")

# Streamlit UI
st.title("PubMed Research Navigator")
st.write("Search PubMed for articles and save the results as an Excel file.")

email = st.text_input("Enter your email (for PubMed access):")
search_term = st.text_input("Enter the general search term:")
mesh_term = st.text_input("Enter an optional MeSH term (leave blank if not needed):")

article_choice = st.selectbox("Select article type:", list(article_types.keys()))
num_articles = st.number_input("Enter the number of articles to fetch:", min_value=1, max_value=1000, value=10)

if st.button("Fetch Articles"):
    if email and search_term:
        query = construct_query(search_term, mesh_term, article_choice)
        articles = fetch_abstracts(query, num_articles, email)
        if articles:
            excel_data = save_to_excel(articles)
            st.download_button(
                label="Download Abstracts",
                data=excel_data,
                file_name="pubmed_articles.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
            year_count = count_articles_by_year(articles)
            if year_count:
                st.write("Publication Year Distribution:")
                plot_publication_years_pie_chart(year_count)
            analyze_disease_frequency(articles)
        else:
            st.write("No articles fetched.")
    else:
        st.write("Please fill in all the required fields.")

st.write("""
    ### Copyright Information
    Copyright (c) 2024 Anurag Pande
""")
