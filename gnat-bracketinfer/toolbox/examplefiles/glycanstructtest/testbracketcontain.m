clc;clear;
a=glycanMLread('TEST1.glycoct_xml');
b=glycanMLread('TEST2.glycoct_xml');
glycanViewer(a);
glycanViewer(b);
ajava=a.glycanjava;
bjava=b.glycanjava;