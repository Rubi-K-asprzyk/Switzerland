TracheoData = Matrice présence absence pour les 1057 plots où j'ai viré les espèces indéterminées. Par contre il reste toujours les aggrégats.
Matrix_LV95_cleaned_V2 = Matrice pres/absence mousse pour les 575 plots
EnvironmentalVariables_...(_Tracheo) = Les matrices environnementales pour les bryo (trachéo). Note que les NA pour les Snow Cover Duration correspondent à une valeur de 50. La colonne EauxStagnante (la 59 ou 47eme suivant la matrice) n'est pas utile car pas de variation
TableSX.. = La description des variables
EnvironmentalVariablesNames = Le nom des prédicteurs + leurs catégories de variables (LC = land Cover, Clim = climatic, Edaphic, Topo = topographic)

Dis-moi si tu as des questions sur le jeu de données. Les coordonnées sont dans le CRS EPSG:2056 = CH1903+/LV95.

- Mnt_Mean renamed into "Z" 

- Initial bryophytes occurence data modified
    - Removal of the agg. sp. sterile ...
        - Bryum_sp
        - Cephaliozella_sp
        - Encalypta_sp
        - Fossombronia_sp
        - Hypnum_sp
        - Leiocolea_sp
        - Marsupella_sp
        - Pellia_sp
        - Pohlia_sp
        - Scapania_sp
        - Schistidium_sterile
        - Weissia sterile
        - Eurhynchiastrum_pulchellum_diversifolium -> Eurhynchiastrum_pulchellum
        - Sciuro_hypnum_X -> Sciurohypnum_X

    - Fusion of Palustriella sulcata and Palustriella falcata (1 seule occurence) into Palustriella sulcata. 