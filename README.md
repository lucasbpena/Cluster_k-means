Authors: 
  - Marcos G. Quiles (@quiles) 
  - Marinalva D. Soares

1º GERAR ESTRUTURAS ALEATÓRIAS
  1. Arquivo .xyz
    - Deixar só as coordenadas atômicas (apagar o número 55, a linha em branco e o elemento)
    - novo arquivo M55
  2. Criar uma pasta com a composição desejada (AnBm)
    - Colar o arquivo M55
    - Colar o script randgenfinal.py
  3. Abrir o randgenfinal.py
    - Mudar o nome do arquivo (M55)
    - Mudar o nome do range (número de estruturas que devem ser geradas)
  4. Rodar o script randgenfinal.py
    - python2 randgenfinal.py
      *adicionar elemento 1 e a composição 1
      *adicionar elemento 2 e a composição 2
2º PROCESSO DE CLUSTERIZAÇÃO
  Na pasta anterior, rodar o script “silscript.py”
    python3 silscript.py 1 AnBm n
  Onde:
    AnBm - pasta onde será rodado o script
    n - número de estruturas que devem ser selecionadas
