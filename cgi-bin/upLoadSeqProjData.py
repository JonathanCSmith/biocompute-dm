#!/usr/bin/python


import cgi

CH = open("../template/cgi_header", "r")
cgiHeader = CH.readlines()
CH.close()

HF = open("../template/template_top.html", "r")
topLines = HF.readlines()
HF.close()

FF = open("../template/template_bottom.html", "r")
bottomLines = FF.readlines()
FF.close()

for ch in cgiHeader:
    print
    ch[:-1]

for he in topLines:
    print
    he[:-1]

print
"<p>Upload csv file:</p>"
print
'<form name="input" enctype="multipart/form-data" action="/cgi-bin/upLoaderSeqProjStage1" method="post">'
print
'<p><input type="file" name="seqDataFile" size="40"><br> <input type="submit" Value="Upload" ></p>'
print
'</form>'
"""
print "<p>or</p>"

print '<form name="pasteinput" action="/cgi-bin/upLoaderSeqProjStage1" method="post">'
print "<p>Paste csv file here:</p>"
print '<textarea name="csv" width="400" height="50" > </textarea>'
print '<p><input type="submit" Value="Upload" ></p>'
print '</form>'
"""
print
"<hr>"
print
"<h2>Important - Sample Name Illegal Characters</h2>"
print
"""<p>A number of characters are not allowed in the sample names(Customer Sample ID) input to casava, these are the space character and the following :<br><b>? ( ) [ ] / \ = + < > : ; " ' , * ^ | & .</b></p>"""

print
"<h2>Upload File Format</h2>"

print
"<p>Files must be CSV format, and the first row should always contain the column headers."

print
"A minimum set of information from each customer project is:"
print
"<ul><li>Lane number(s)</li><li>Customer Sample ID</li><li>Index sequence(s) (if multiplexed)</li></ul>"
print
"<p>The following is an example of the layout of the csv file:</p>"
print
"<table border='1' style='background-color:lightblue;' >"
print
"<tr><th>ProjectID</th><th>Sequencing Concentration</th><th>PhiXspiked</th><th>Lane</th><th>Customer Sample ID</th><th>Index sequence</th></tr>"
print
"<tr><td>P303</td><td>10pM</td><td>Y</td><td>1</td><td>BR1</td><td>TAAGGCGA</td></tr>"
print
"<tr><td>P303</td><td>10pM</td><td>Y</td><td>1</td><td>BR3</td><td>CGTACTAG</td></tr>"
print
"<tr><td>P303</td><td>10pM</td><td>Y</td><td>1</td><td>BR15_A</td><td>AGGCAGAA</td></tr>"
print
"<tr><td>P303</td><td>10pM</td><td>Y</td><td>1</td><td>CON_1</td><td>TCCTGAGC</td></tr>"
print
"<tr><td>P303</td><td>10pM</td><td>Y</td><td>1</td><td>BR4</td><td>GGACTCCT</td></tr>"
print
"<tr><td>P303</td><td>10pM</td><td>Y</td><td>1</td><td>BR9</td><td>TAGGCATG</td></tr>"
print
"<tr><td>P303</td><td>10pM</td><td>Y</td><td>1</td><td>BR17</td><td>CTCTCTAC</td></tr>"
print
"<tr><td>P303</td><td>10pM</td><td>Y</td><td>1</td><td>BR18</td><td>CAGAGAGG</td></tr>"
print
"<tr><td>P303</td><td>10pM</td><td>Y</td><td>1</td><td>BR19</td><td>GCTACGCT</td></tr>"
print
"<tr><td>P303</td><td>10pM</td><td>Y</td><td>1</td><td>CON_2</td><td>CGAGGCTG</td></tr>"
print
"<tr><td>P303</td><td>10pM</td><td>Y</td><td>1</td><td>BR20</td><td>AAGAGGCA</td></tr>"
print
"<tr><td>P303</td><td>10pM</td><td>Y</td><td>1</td><td>BR21</td><td>GTAGAGGA</td></tr>"
print
"<tr><td>P303</td><td>14pM</td><td>Y</td><td>2</td><td>CC1</td><td>TAAGGCGA</td></tr>"
print
"<tr><td>P303</td><td>14pM</td><td>Y</td><td>2</td><td>CC3</td><td>CGTACTAG</td></tr>"
print
"<tr><td>P303</td><td>14pM</td><td>Y</td><td>2</td><td>CC15_B</td><td>AGGCAGAA</td></tr>"
print
"<tr><td>P303</td><td>14pM</td><td>Y</td><td>2</td><td>CON_3</td><td>TCCTGAGC</td></tr>"
print
"<tr><td>P303</td><td>14pM</td><td>Y</td><td>2</td><td>CC4</td><td>GGACTCCT</td></tr>"
print
"<tr><td>P303</td><td>14pM</td><td>Y</td><td>2</td><td>CC9</td><td>TAGGCATG</td></tr>"
print
"<tr><td>P303</td><td>14pM</td><td>Y</td><td>2</td><td>CC17</td><td>CTCTCTAC</td></tr>"
print
"<tr><td>P303</td><td>14pM</td><td>Y</td><td>2</td><td>CC18</td><td>CAGAGAGG</td></tr>"
print
"<tr><td>P303</td><td>14pM</td><td>Y</td><td>2</td><td>CC19</td><td>GCTACGCT</td></tr>"
print
"<tr><td>P303</td><td>14pM</td><td>Y</td><td>2</td><td>CON_4</td><td>CGAGGCTG</td></tr>"
print
"<tr><td>P303</td><td>14pM</td><td>Y</td><td>2</td><td>CC20</td><td>AAGAGGCA</td></tr>"
print
"<tr><td>P303</td><td>14pM</td><td>Y</td><td>2</td><td>CC21</td><td>GTAGAGGA</td></tr>"
print
"</table>"

for fo in bottomLines:
    print
    fo[:-1]
