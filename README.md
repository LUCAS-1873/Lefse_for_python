# python实现Lefse的脚本
By: 中南大学 人体与微生物学 赵一鸣   


这是一个可以通过python3版本实现Lefse的python脚本。可以在linux下运行。

Pandas == 1.1.3  
Numpy == 1.19.2  
Scipy == 1.5.2  
Sklearn==0.23.2  

  
**用法：**  

usage: lefse_by_Lucas.py [-h] -s SPTAB -g GROUP [-l LDA_TH] -o OUTPUT  

optional arguments:  
  -h, --help            show this help message and exit  
  -s SPTAB, --sptab SPTAB  
                        sptab.txt  
  -g GROUP, --group GROUP  
                        group.txt  
  -l LDA_TH, --lda_th LDA_TH  
                        default=2.0  
  -o OUTPUT, --output OUTPUT  
                        out_dir  
  
python lefse_by_Lucas.py -s sptab.txt -g group.txt  -o out_dir  
  
  
注意：  
1. 输入文件格式为CSV格式；请参考sample_data文件夹中的参考文件格式。  
2. 若返回结果为空，请调小LDA的阈值再尝试。  
3. 本脚本输入的数据并未进行任何归一化或标准化。如需要请自行调整输入数据格式。  
4. 本脚本只考虑一层级关系，即只考虑一分组关系，若该分组下还有层级则无法完成。  
5. 引用请告知本人谢谢！  
   
特别感谢中科院地球化学所吴求生在本脚本中做出的调试及贡献！  

若需要使用 或有任何疑问 请联系 yimingzhao970904@163.com 
