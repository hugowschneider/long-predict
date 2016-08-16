# long-predict

## Instalation
```
wget https://github.com/hugowschneider/long-predict/archive/v0.9.tar.gz
tar zxvf v0.9.tar.gz
cd long-predict-0.9
cmake .
make
```
## Usage
```
Usage: ./long_predict -c <model config file> -i <fasta file> [-d <both|-|+> ] [-s <size>]
	-c	Model config file. File specifying the libSVM model and the used attributes
	-i	Input Fasta file for prediction. This file should be a plain text or a gzip fasta
		file
	-d	Direction. The direction of the fasta file sequences for prediction. '+' will
		read the sequences as is. '-' will use the complementary sequence. And both will
		use both options. This option is used for all sequences in the file. Default '+'
	-s	Size limit. This attribute ignore sequences shorter than this limit. Default 200

Model Config File is a plain text file containing the following attributes:
	modelFile	The path to the model file. It can be relative to the config file or
				a absolute path
	attributes	The list of attributes used in the model training. This attributes are valid
				nucleotide frequencies, for example 'aa' and 'atc', and the values 'ol' for
				first ORF lenght and 'op' for first ORF percentage of the corresponding transcript
				length
Ex.:
modelFile=human.model
attributes=aa,aaa,ac,aca,acg,op

```

To classify all sequences in the file `human.fa.gz` using the configuration file located in 
`models/human.model.conf`: 

```
./long_predict -c models/human.model.conf -i human.fa.gz
```

This outputs the following:

```
Predicted lncRNAs in complementary sequence: 40.00% (40/100)
```

And it creates the csv file `human.fa.gz.csv` containing all probalities for the prediction:
```
ID,Size,Classification,Probability lncRNA,Probability PCT
"ENST00000437894",227,lncRNA,0.9973536549711840,0.0026463450288160
...
```


## Copyright
The work herein is Copyright 20013--2016 Hugo Wruck Schneider and Universidade de Brasília (UnB). **No rights are given to reproduce or modify this work**.

This work uses libSVM for model training and prediction:

Chih-Chung Chang and Chih-Jen Lin, LIBSVM : a library for support vector machines. ACM Transactions on Intelligent Systems and Technology, 2:27:1--27:27, 2011. Software available at [http://www.csie.ntu.edu.tw/~cjlin/libsvm](http://www.csie.ntu.edu.tw/~cjlin/libsvm)