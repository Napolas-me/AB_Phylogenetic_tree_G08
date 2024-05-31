# AB_Phylogenetic_tree_G08

Para a implementação em Python é necessario instalar a biblioteca Biopython. Esta instalação pode ser realizada através de um comando de instalação de dependencias para python como:
```
  pip install biopython
```
O Stript 'concat_fast_files.py' não necessita de qualquer atualização de nomes de ficheiros, apenas será necessario correr o programa atraves de um terminal ou IDE dentro da pasta "python" do repositorio.
Este script irá gerar um ficheiro, concatenated.fasta dentro da pasta "Species", com a concatenação de todos os ficheiras FASTA presentes dentro da pasta "Species/data". Se for do interresse utilizar um dataset diferente ou maior basta modificar o conteudo da pasta "Species/data". 

Dentro do script de python "read_msa_to_tree" é necessário atualizar o nome do ficheiro a ser utilizado para a geração das árvores. Se necessario os nomes das imagens criadas também podem ser modificados no código.
