# Structure générale du code AnalyzeJetsGrid.cpp

## 📋 Vue d'ensemble

Ce code C++ fait partie du projet **ALICE FOCAL** (Forward Calorimeter for ALICE experiment au CERN). Il analyse les jets de particules dans une grille d'événements, effectue l'appariement des jets au niveau détecteurs avec les jets au niveau particules issues de la simulation Monte-Carlo, et stocke les résultats dans un fichier ROOT pour une analyse ultérieure.

---

## 🔧 Dépendances et bibliothèques

### Framework ROOT (analyse de données)
- `TFile.h` - Gestion des fichiers ROOT
- `TTree.h` - Arbres de données
- `TClonesArray.h`, `TObjArray.h` - Conteneurs de données
- `TH1.h`, `TH2.h`, `TH3.h` - Histogrammes
- `TProfile.h` - Profils d'histogrammes
- `TParticle.h` - Représentation des particules
- `Riostream.h` - Flux d'entrée/sortie

### Bibliothèques ALICE (simulation et détection)
- `AliRunLoader.h` - Gestion des événements et de la pile MC
- `AliStack.h` - Pile de particules générées
- `AliFOCALCluster.h` - Clusters du détecteur FOCAL

### FastJet (reconstruction de jets)
- `fastjet/PseudoJet.hh` - Classe représentant les jets
- `fastjet/ClusterSequenceArea.hh` - Clustering avec calcul d'aire
- `fastjet/GhostedAreaSpec.hh` - Spécification de fantômes pour le calcul d'aire

---

## 📦 Déclarations de fonctions (forward declarations)

```cpp
Float_t eta(TParticle *);           // Calcule la pseudorapidité η
Float_t phi(TParticle *);           // Calcule l'angle azimutal φ
void GetMomentum(...);              // Calcule la 4-quantité de mouvement
void GetGeometricalMatchingLevel(...); // Distance géométrique entre jets
void DoJetMatching(...);            // Appariement des jets
```

---

## 🏗️ Classes personnalisées

### 1. **Classe `JetMatchingParams`**
Hérite de `fastjet::PseudoJet::UserInfoBase`

**Objectif:** Stocker les informations d'appariement des jets (indices et distances)

**Membres privés:**
- `mIndex1` - Indice du jet apparié le plus proche
- `mIndex2` - Indice du jet apparié deuxième plus proche
- `mDistance1` - Distance au jet le plus proche
- `mDistance2` - Distance au deuxième jet le plus proche

**Méthodes publiques:**
- `index1()`, `index2()` - Accès aux indices
- `distance1()`, `distance2()` - Accès aux distances
- `setJetIndexDistance1()`, `setJetIndexDistance2()` - Définition des indices et distances

---

### 2. **Classe `SW_JetArea`**
Hérite de `fastjet::SelectorWorker`

**Objectif:** Sélecteur personnalisé pour filtrer les jets basé sur leur aire

**Contexte:** L'aire représente la surface occupée dans l'espace de phase (η, φ) par le jet. Elle est utilisée pour:
- La soustraction du bruit de fond (background subtraction): `p_T^sub = p_T - ρ × Area`
- Le jet grooming (élagage)

**Méthodes:**
- `SW_JetArea(double minArea)` - Constructeur avec aire minimale
- `pass(const fastjet::PseudoJet &p)` - Retourne vrai si l'aire du jet > minArea
- `description()` - Description du sélecteur

**Fonction associée:**
```cpp
fastjet::Selector SelectorMinArea(const double &minArea)
```
Crée un sélecteur pour une aire minimale donnée.

---

## 🎯 Fonction principale: `AnalyzeJetsGrid()`

Cette fonction est le cœur du programme. Voici son flux général:

### **Étape 1: Initialisation des paramètres**
```
- Ecut = 2.0 GeV (seuil d'énergie)
- Masses: π⁰ (0.135 GeV), π⁺ (0.139 GeV)
- Rayons de cône: R = 0.2, 0.4, 0.6
```

### **Étape 2: Ouverture des fichiers**
- Lecture des clusters du fichier `focalClusters.root`
- Création du fichier de sortie `analysisJets.root`

### **Étape 3: Création des arbres ROOT**
Deux arbres de données sont créés:

#### **`jetTree`** (jets reconstruits)
Branches:
- `ievt` - Numéro d'événement
- `jetE`, `jetpT` - Énergie et impulsion transverse
- `jetpT_sub` - Impulsion transverse après soustraction du bruit
- `jetPhi`, `jetEta`, `jetRap` - Coordonnées (φ, η, rapidité)
- `jetE_match`, `jetpT_match` - Propriétés du jet apparié (vérité MC)
- `jetPhi_match`, `jetEta_match`, `jetRap_match` - Coordonnées du jet apparié
- `jet_distmatch` - Distance d'appariement
- `jetR` - Rayon du cône utilisé

#### **`truthjetTree`** (jets au niveau Monte-Carlo)
Branches:
- `truthjetE`, `truthjetpT` - Énergie et impulsion
- `truthjetPhi`, `truthjetEta`, `truthjetRap` - Coordonnées
- `truthjetR` - Rayon du cône

### **Étape 4: Chargement du générateur d'événements (AliRunLoader)**
- Charge la pile Monte-Carlo (stack)
- Charge les informations cinématiques
- Récupère les primaires physiques

### **Étape 5: Définitions de jets**
Crée des objets `JetDefinition` pour:
- Algorithme **anti-kt** (reconstruction): utilisé pour les jets reconstruits
- Algorithme **kt** (clustering hiérarchique): utilisé pour le calcul du bruit

Configuration avec aire (ghost area spec pour background subtraction):
```cpp
Double_t ghost_maxrap = 6.0
AreaDefinition area_def(fastjet::active_area, area_spec)
```

### **Étape 6: Boucle sur les événements**

Pour chaque événement:

#### **6a. Récupération de la pile Monte-Carlo**
```
- Itère sur toutes les particules générées
- Récupère la position du vertex (de la 6ème particule)
- Sélectionne les primaires physiques
- Garde uniquement les pseudo-rapidité tel que 3.2 < η < 5.8 (région forward FOCAL)
- Crée un vecteur PseudoJet pour FastJet
```

#### **6b. Récupération des clusters du détecteur**
```
- Charge les clusters FOCAL du fichier externe
- Combine les clusters ECAL (électromagnétique)
- Ajoute les clusters HCAL (hadronique) du segment 6
- Rejette les clusters avec énergie anormale (> 6500)
- Crée un vecteur PseudoJet pour FastJet
```

#### **6c. Reconstruction des jets (pour chaque rayon R)**

**Pour chaque rayon de cône (R = 0.2, 0.4, 0.6):**

1. **Clustering des jets reconstruits** (à partir des clusters du détecteur)
   ```cpp
   ClusterSequenceArea clust_seq(input_Clusters, jet_def[iR], area_def)
   inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin=1.0 GeV))
   ```

2. **Clustering des jets vérité** (à partir des particules MC)
   ```cpp
   ClusterSequenceArea clust_seq_truth(input_Particles, jet_def[iR], area_def)
   inclusive_jets_truth = sorted_by_pt(clust_seq_truth.inclusive_jets(ptmin=1.0 GeV))
   ```

3. **Calcul du bruit de fond (ρ)**
   - Utilise un cône perpendiculaire pour éviter le jet signal
   - Position: φ = φ(jet₀) + π/2, η = η(jet₀)
   - `Rho = PerpConePt / (π × R²)`

4. **Initialisation des structures de matching**
   ```cpp
   - Chaque jet reçoit un index utilisateur
   - Chaque jet est associé à un objet JetMatchingParams vide
   ```

5. **Appariement des jets** (voir fonction `DoJetMatching()`)

6. **Remplissage des arbres ROOT**
   ```
   - Pour chaque jet reconstruit:
     * Remplir les propriétés du jet
     * Remplir les propriétés du jet appié (vérité MC)
     * Remplir la distance d'appariement
   - Pour chaque jet vérité:
     * Remplir les propriétés du jet MC
   ```

### **Étape 7: Fermeture et sauvegarde**
```cpp
- Écriture des arbres dans le fichier ROOT
- Fermeture des fichiers
- Nettoyage des structures (Delete)
```

---

## 🔗 Fonctions utilitaires

### 1. **`Float_t eta(TParticle *part)`**
Calcule la **pseudorapidité η** d'une particule
```
η = -ln(tan(θ/2)) où θ est l'angle polaire
```
**Entrée:** Pointeur sur particule TParticle  
**Sortie:** Pseudorapidité (retourne 9999 si pt = 0)

---

### 2. **`Float_t phi(TParticle *part)`**
Calcule l'**angle azimutal φ** d'une particule
```
φ = atan2(Py, -Px)
```
**Entrée:** Pointeur sur particule TParticle  
**Sortie:** Angle azimutal en radians

---

### 3. **`void GetMomentum(TLorentzVector &mom, const Float_t *vertex, AliFOCALCluster *foClust, Float_t mass)`**

Calcule la **4-quantité de mouvement** d'un cluster

**Paramètres:**
- `mom` - TLorentzVector qui recevra le 4-vecteur
- `vertex` - Position du vertex primaire
- `foClust` - Cluster du détecteur FOCAL
- `mass` - Masse de la particule (par défaut 0)

**Algorithme:**
```
1. Récupère l'énergie E du cluster
2. Calcule l'impulsion: p = √(E² - m²)
3. Récupère la position du cluster (x, y, z)
4. Soustrait le vertex pour obtenir la direction
5. Normalise et multiplie par p:
   Px = p × (x-vx)/r
   Py = p × (y-vy)/r
   Pz = p × (z-vz)/r
   où r = √((x-vx)² + (y-vy)² + (z-vz)²)
```

---

### 4. **`void DoJetMatching(std::vector<fastjet::PseudoJet> &jetArray1, std::vector<fastjet::PseudoJet> &jetArray2)`**

Effectue l'**appariement bidirectionnel des jets** entre deux collections

**Paramètres:**
- `jetArray1` - Jets reconstruits (sera modifiée)
- `jetArray2` - Jets vérité MC (sera modifiée)

**Algorithme:**

#### **Phase 1: Calcul des distances**
```
Pour chaque jet1 dans jetArray1:
  Pour chaque jet2 dans jetArray2:
    Calcule distance = Δ R = √((Δη)² + (Δφ)²)
    Met à jour les deux meilleurs appariements pour jet1
    Met à jour les deux meilleurs appariements pour jet2
```

#### **Phase 2: Réconciliation des appariements mutuels**
```
Pour chaque jet1:
    Si son jet1 le plus proche (jet2) ne le pointe pas en retour:
    Essaie de l'assigner à son 2ème choix
    Sinon, marque comme sans appariement
```

**Résultat:** Chaque jet reçoit l'indice et la distance de son jet apparaié le plus proche

---

### 5. **`void GetGeometricalMatchingLevel(fastjet::PseudoJet& jet1, fastjet::PseudoJet& jet2, Double_t &d)`**

Calcule la **distance géométrique en (η, φ)** entre deux jets

```
Δη = η₂ - η₁
Δφ = φ₂ - φ₁  (normalisé à [-π, π])
d = √((Δη)² + (Δφ)²)
```

**Note:** Cette fonction n'est pas utilisée par défaut; elle est commentée dans `DoJetMatching()` au profit de `delta_R()` (qui utilise la rapidité)

---

## 📊 Flux de données général

```
┌─────────────────────────────────────┐
│ Fichier: root_archive.zip#galice.root
│ (Événements MC générés)
└──────────────┬──────────────────────┘
               │
               ▼
     ┌────────────────────┐
     │ AliStack           │
     │ (Particules primaires)
     └────────────────────┘
               │
               ▼
     ┌────────────────────┐
     │ FastJet Clustering │
     │ (anti-kt R=0.2,0.4,0.6)
     └────────────────────┘
               │
               ▼
     ┌────────────────────┐
     │ inclusive_jets_truth
     │ (Jets vérité MC)
     └────────────────────┘
               │
               ▼
┌──────────────────────────────────────┐
│ Fichier: focalClusters.root
│ (Clusters détecteur reconstruits)
└──────────────┬───────────────────────┘
               │
               ▼
     ┌────────────────────┐
     │ FastJet Clustering │
     │ (anti-kt + area)
     └────────────────────┘
               │
               ▼
     ┌────────────────────┐
     │ inclusive_jets
     │ (Jets reconstruits)
     └────────────────────┘
               │
               ▼
     ┌────────────────────┐
     │ DoJetMatching      │
     │ (Appariement)
     └────────────────────┘
               │
               ▼
     ┌────────────────────┐
     │ Calcul du bruit ρ  │
     │ (cône perpendiculaire)
     └────────────────────┘
               │
               ▼
     ┌────────────────────┐
     │ Remplissage arbres │
     │ ROOT
     └────────────────────┘
               │
               ▼
     ┌────────────────────┐
     │ analysisJets.root
     │ (Résultats)
     └────────────────────┘
```

---

## 🎓 Concepts clés

| Concept | Explication |
|---------|------------|
| **Jet** | Collection colimée de particules (ou clusters) reconstruits via algorithme de clustering |
| **Anti-kt** | Algorithme de clustering qui privilégie les grandes structures (jets) |
| **Pseudorapidité (η)** | Variable cinématique invariante par boost longitudinal: η = -ln(tan(θ/2)) |
| **Rapidité (rap)** | Variable cinématique complètement invariante: rap = ½ ln((E+Pz)/(E-Pz)) |
| **ΔR** | Distance dans l'espace (η, φ): ΔR = √((Δη)² + (Δφ)²) |
| **Aire du jet** | Surface occupée dans (η, φ), utilisée pour background subtraction |
| **Background subtraction** | Soustraction du bruit: pT_sub = pT - ρ × Area |
| **Vérité MC** | Particules/jets au niveau simulation Monte-Carlo (avant détecteur) |
| **Clusters reconstruits** | Signaux mesurés dans le détecteur FOCAL |

---

## 📌 Fichiers d'entrée/sortie

### **Entrée**
- `root_archive.zip#galice.root` - Événements MC (pile de particules)
- `focalClusters.root` - Clusters reconstruits du détecteur FOCAL

### **Sortie**
- `analysisJets.root` - Arbre des jets reconstruits et appariés

---

## 🚀 Résumé du processus

1. **Charger** les événements MC et les clusters détecteur
2. **Reconstruire** les jets vérité à partir des particules MC
3. **Reconstruire** les jets à partir des clusters détecteur
4. **Apparier** automatiquement les jets reconstruits avec les jets vérité
5. **Calculer** le bruit de fond pour la soustraction
6. **Stocker** tous les résultats dans un arbre ROOT pour analyse ultérieure