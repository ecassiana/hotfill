#! /bin/bash

# Runs {LaTex} on an internal document in order to force the
# creation of PK files for small versions of the default fonts, 
# so that hacks.plot_text can work with small {fsize} parameter.

jobname="$$"
texfile="${jobname}.tex"
dvifile="${jobname}.dvi"
pdffile="${jobname}.pdf"
errfile="${jobname}.err"
logfile="${jobname}.log"

pushd /tmp/

cat <<EOF > ${texfile}
\nonstopmode
\documentclass{article}
\usepackage{lmodern}
\usepackage{anyfontsize}

\newcommand{\trysize}[1]{%
\hrule
\fontsize{#1}{#1}\selectfont
  \begin{center}
    foo 

    \textbf{bf:foo} 
    \textit{it:foo}
    \textsf{sf:foo} 
    \textsl{sl:foo}
    \texttt{tt:foo}

    \textit{\textbf{it:bf:foo}}
    \textsf{\textbf{sf:bf:foo}} 
    \textsl{\textbf{sl:bf:foo}}
    \texttt{\textbf{tt:bf:foo}}

    \textbf{\textit{bf:it:foo}} 
    \textsf{\textit{sf:it:foo}} 
    \textsl{\textit{sl:it:foo}}
    \texttt{\textit{tt:it:foo}}

    \textsf{\textsl{sf:sl:foo}} 
    \textbf{\textsl{bf:sl:foo}} 
    \texttt{\textsl{tt:sl:foo}}
    \textit{\textsl{it:sl:foo}}

    \textsf{\textbf{\textit{sf:bf:it:foo}}} 
    \texttt{\textbf{\textit{tt:bf:it:foo}}}
    \texttt{\textbf{\textsl{tt:bf:sl:foo}}}
  \end{center}
\hrule
}

\begin{document}
\trysize{1}
\trysize{5}
\trysize{37}
\end{document}
EOF

which latex

latex ${jobname} > ${errfile}
texfailed=$?
echo "texfailed = ${texfailed}"
if [[ ( ${texfailed} -eq 0 ) && ( -s ${dvifile} ) ]]; then
  dvipdf  ${dvifile} ${pdffile}
  evince ${pdffile}
fi
echo "================================================================================" 1>&2
cat ${errfile} 1>&2
echo "================================================================================" 1>&2
cat ${logfile} 1>&2
echo "================================================================================" 1>&2

popd
