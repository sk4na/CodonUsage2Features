# CodonUsage2Features
Bioinformatic workflow used to elucidate the possible relationship between codon usage data, as a proxy of mRNA translational speed,
and different protein features mapped to the protein sequences by UniProt, in multiple organisms. Code developed as a part of my masters thesis.

######################## Uso de los scripts ########################
## Nota: en la carpeta 'data', faltan los archivos de uniprot ######
## debido a que algunos superan el tamaño permitido por GitHub. ####
## Habría que incorporarlos manualmente en caso de querer ejecutar #
## el flujo de trabajo y replicar los resultados obtenidos #########
####################################################################

#### workflow_launch.sh ####
Este es el script maestro, un script de shell utilizado para ejecutar el flujo de trabajo en su totalidad. No acepta argumentos, simplemente se ejecuta
en la linea de comandos y el flujo de trabajo se lleva a cabo. Por dentro, este script ejecuta los siguientes scripts:
  
  ## seq2features.sh ##
  Script de bash con dos posibles modos de uso:
    opción -c: modo 'create', toma como argumentos el fichero de KEGG, el de UniProt, el de codon usage y el nombre del organismo para el cual se
    va a hacer el análisis. El output es el fichero con las proteínas filtradas y anotadas a nivel de cada codón que se explica en la metodología
    del TFM. En este modo, el script que se ejecuta por dentro es el script de python 'get_annot_sequences.py'.
    opción -v: modo 'visualize', el cual ejecuta el script de python 'visualize_proteins.py' utilizado para la visualización de proteínas concretas
    generando gráficas en formato SVG como las mostradas en la sección 4.4 de los resultados. Este script toma como argumentos:
    - el archivo generado con la opción '-c' de seq2features.sh
    - el archivo de uniprot
    - una lista separada por comas de AC de UniProt para las proteínas que se quieran visualizar (argumento -p)
    - la ventana de posiciones que se desee emplear para el ajuste de los datos de RSCU mostrados (-w).
    
  ## analyze_features.py ##
  Script de python encargado de crear los perfiles medios de valores RSCU y de datos sobre la posición de codones óptimos/raros (RPCU),
  llevar a cabo el análisis estadístico descrito en la metodología, y representar gráficamente los perfiles medios obtenidos con respecto a la característica
  proteica analizada. Además, escribe en ficheros todos los perfiles de datos de codon usage creados para cada proteína del organismo en ficheros
  ordenados por nombres y en carpetas ordenadas. El script toma como argumentos: 
  - el archivo generado por la opción -c de 'seq2features.sh'
  - un output path, que por defecto es el path de ejecución del script
  - un fichero con una feature de UniProt por fila, tantas como se quieran analizar (-f)
  - una string con el tipo de datos que se va a analizar, 'RSCU' o 'codon_type' (-d)
  - el número de posiciones para la ventana de ajuste de los datos a representar en la gráfica (-w) y el nombre del organismo (-o).
  
  ## perform_all_organism_ttest.py ##
  Script de python utilizado para el análisis conjunto de los perfiles obtenidos de todos los organismos estudiados. Funcionamiento
  y argumentos muy similares a 'analyze_features.py', solo que en este caso toma 2 archivos de input:
  - archivo que contiene todos los perfiles obtenidos de todos los organismos para una feature concreta, generado por el script 'workflow_launch.sh'
  - archivo que contiene todos los perfiles aleatorios creados para todos los organismos, generado por el script 'workflow_launch.sh'.
  
  ## all_organism_graphs.py ##
  Script de python utilizado para generar las figuras 3 y 8, mostradas en los resultados. El script toma como argumentos:
  - Archivo que contiene todos los perfiles obtenidos de todos los organismos para una feature concreta
  - un output path, que por defecto es el path de ejecución del script
  - una string con el tipo de datos que se va a analizar, 'RSCU' o 'codon_type' (-d)
  - el número de posiciones para la ventana de ajuste de los datos a representar en la gráfica (-w) y el nombre del organismo (-o).
  
  ## compare_data_types.py ##
  Script de python utilizado para generar la figura 13. Toma como argumentos:
  - el archivo generado por la opción -c de 'seq2features.sh'
  - archivo con los perfiles de valores RSCU de E.coli para la feature que se quiera analizar
  - archivo con los perfiles de codon_type de E.coli para la feature que se quiera analizar
  - archivo con los datos de [tRNA] de E.coli, cuya obtención se menciona en la metodología
  - un output path, que por defecto es el path de ejecución del script
  - string con el nombre de la feature que se va a analizar (-f)
  - el número de posiciones para la ventana de ajuste de los datos a representar en la gráfica (-w) y el nombre del organismo (-o).


Adicionalmente, en el repositorio se ha añadido un archivo zip que contiene todas las gráficas de perfiles medios de datos de codon usage
generadas durante el desarrollo del TFM.
  
