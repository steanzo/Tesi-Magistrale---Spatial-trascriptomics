# Spatial-Trascriptomics
Nel mio progetto di tesi ho sviluppato un modello statistico avanzato per analizzare dati di trascrittomica spaziale, una tecnologia che permette di studiare l’espressione genica preservando la posizione delle cellule all’interno del tessuto. Come riportato nel documento, “l’obiettivo di questo lavoro è sviluppare un modello statistico che affronti efficacemente sovradispersione e dipendenza spaziale nei dati omici”.

Per farlo ho integrato:
- Distribuzione Binomiale Negativa per gestire la sovradispersione tipica dei dati RNA-seq.
- Campo Gaussiano per modellare la correlazione tra spot vicini.
- Approccio SPDE (Stochastic Partial Differential Equations) per stimare in modo efficiente il campo spaziale.
- Inferenza Bayesiana tramite INLA, che ha permesso stime rapide e accurate anche su dataset complessi.

Il modello è stato applicato a dati Visium (10x Genomics) provenienti da tessuto muscolare murino, integrati con informazioni da immagini di immunofluorescenza. L’obiettivo biologico era valutare l’effetto di un trattamento su geni coinvolti nella degradazione muscolare, come Trim63, Foxo1 e Foxo3. Come riportato nel testo, “il modello ha evidenziato effetti significativi del trattamento su 196 geni, principalmente legati a processi immunitari”, mentre i tre geni di interesse non hanno mostrato l’effetto atteso.

Il lavoro ha permesso di:
- costruire un pipeline completa di preprocessing, normalizzazione e integrazione immagine–RNA;
- confrontare modelli spaziali alternativi (BN, Poisson, GMRF, SPDE);
- evidenziare pattern spaziali di espressione e aree di sovra/sotto-espressione nei tessuti;
- proporre miglioramenti futuri, come l’uso di prior spike-and-slab per aumentare la potenza nei test multipli.

## Descrizione
Questo repository raccoglie esclusivamente i codici utilizzati nel mio lavoro di tesi: **"Analisi di dati omici con dipendenza spaziale: un modello binomiale negativo con campo gaussiano"**, svolto presso l'Università degli Studi di Padova, per il corso magistrale in Scienze Statistiche. L'obiettivo è organizzare e documentare in maniera ordinata tutto il codice implementato.

## Struttura del Progetto
- **allineamento**: Contiene tutti gli script utilizzati per l'allineamento di immagini
- **estrazione**: Contiene i codicie relativi all'estrazione d'informazione dalle immagini
- **adattamento modello**: Contiene i codici utilizzati per l'adattamento del modello
- **codici utili** vengono anche presentati dei codici per l'estrazione di grafici dal modello INLA

## Programmi e librerie utilizzate
Per l'allineamento e l'estrazione del'informazione dai fail `zarr` contenenti l'immagini dei tassuti con immunofluorescenza e i livelli di espressione genica è stato utilizzato il linguaggio di programmazione **python** e i seguenti pacchetti: **spatialdata**, **napari** e **squidpy**.
Per l'adattamento del modello è stato utilizzato il programma **R** e il pacchetto principale è **INLA**.



