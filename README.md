# 🧬 MetaBee

**MetaBee** – A reproducible workflow for genome-scale metabolic reconstruction and interaction analysis in the *Apis mellifera* gut microbiome.

---

## 🐝 Overview

**MetaBee** is a computational workflow designed to reconstruct, curate, and analyze **genome-scale metabolic models (GEMs)** for the honeybee (*Apis mellifera*) gut microbiota.

Developed as part of a Master's project conducted in the **Hatzimanikatis Lab** at **EPFL** and in collaboration with the **Engel Lab** at the **University of Lausanne**, the pipeline integrates existing genome annotation and model reconstruction tools into a unified, reproducible process. It enables researchers to build high-quality GEMs and explore potential metabolic interactions within synthetic or natural bee gut microbial communities.

While originally tailored to a **25-member synthetic bee gut community**, MetaBee can be adapted to other microbiomes.

---

## ⚙️ Key Features

- 🧩 **Automated GEM reconstruction** using tools like RAVEN 2.0 
- 🔍 **Integrated genome annotation** via Bakta
- 🧠 **Semi-automated curation** for improving model completeness and consistency using the NICEgame workflow 
- 🌐 **Interaction analysis** to infer potential metabolic exchanges and dependencies
- 🧪 **Reproducible workflow** implemented in Python, MATLAB and integrated into a custom Snakemake pipeline
- 📊 **Modular design** – easily extend or replace components for different datasets

---

## 🧰 Technologies & Dependencies

| Category | Tools |
|-----------|--------|
| **Annotation** | Bakta |
| **Model reconstruction** | RAVEN 2.0 |
| **Modeling frameworks** | COBRApy, nicegamepy, pyTFA |
| **Workflow management** | Snakemake |
| **Utilities** | Python ≥3.9, pandas, numpy, Jupyter, ... |

---

## 🚀 Installation

1. **Clone the repository**

   ```bash
   git clone https://github.com/douhan-wicht/meta-bee.git
   cd meta-bee
2. **TBD**
