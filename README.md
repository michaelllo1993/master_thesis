Automatyzacja status:
- 01_ensembl_peptide_seq_cut.pl (OK. Na wejsciu organizm, wystarczy petla wywolujaca skrypt dla kazdego zadanaego oganizmu)
- 02_run_signalp.pl (OK. Na wejsciu organizm i ilosc plikow wygenerowanych przez 01... W przyszlosci można połączyć z 01...)
- 03_SignalP_positives.pl (OK. Na wejściu organizm, również można połączyć z poprzednimi)
- 04_SigP4_analyzer.pl (OK. Na wejsciu organizm, również można połączyć z poprzednimi)
- shell command: cat sigpL_ALL*/sigp_ALL>>sigpL_ALL.out/sigp_ALL.out.(OK. Wystarczy dodać do pipelinu)
- 05_Needle_wrappers (WYMAGA ROBOTY. Niepotrzebnie jest robiony jest alignmnet człowieka ze wszytskimi organizmami na raz. Zmienić na alignment dla pary podanej na wejściu)

