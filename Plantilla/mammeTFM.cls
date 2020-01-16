%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mammeTFM.cls
%
% Style for the master thesis in MAMME
% This file is a modification of the style file of the 
% journal reports@scm
% to facillitate the publication of the master thesis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\LoadClass[11pt,twoside]{article}
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{reportsascm}[2014/06/25 Estil de la revista reportsascm]
\RequirePackage{xcolor}
\RequirePackage{graphicx}
\RequirePackage{hyperref}
\hypersetup{citecolor=blue}
\RequirePackage{lastpage}
\RequirePackage{sfmath}
\definecolor{blau}{RGB}{43,129,173}
\definecolor{gris}{RGB}{100,100,100}
%Marges
\usepackage[a4paper,left=20mm,top=18mm,right=20mm,bottom=18mm,includeheadfoot,headheight=24mm, headsep=0mm, footskip=15mm]{geometry}
\makeatletter
%Parametres de la pagina
\hypersetup{colorlinks=true,urlcolor=black}
%Definim les variables numero, any, pagin, pagfin, doi, ...
\def\l@any{????}
\def\l@numero{??}
\def\l@pagin{1}
\def\l@pagfin{\pageref*{LastPage}}
\def\l@doi{???????????}
\def\l@abstractcat{Sense abstract}
\def\@nymsc{2010}
\def\@rebut{\today}
\def\@acceptat{\today}
\newcommand{\any}[1]{\def\l@any{#1}}
\newcommand{\numero}[1]{\def\l@numero{#1}}
\newcommand{\pagin}[1]{\def\l@pagin{#1}\setcounter{page}{#1}}
\newcommand{\pagfin}[1]{\def\l@pagfin{#1}}
\newcommand{\doi}[1]{\def\l@doi{#1}}

\def\@msc{}
\def\@graiments{}
\def\@paraulesclau{}
\def\l@abstractceng{No abstract}
\def\l@authorStudent{No abstract}
\def\l@supervisors{No abstract}
\newcommand{\abstracteng}[1]{\def\l@abstracteng{#1}}
\newcommand{\authorStudent}[1]{\def\l@authorStudent{#1}}
\newcommand{\supervisors}[1]{\def\l@supervisors{#1}}
\newcommand{\monthYear}[1]{\def\l@monthYear{#1}}
\newcommand{\agraiments}[1]{\def\@graiments{#1}}
\newcommand{\paraulesclau}[1]{\def\@paraulesclau{#1}}
\newcommand{\titol}[1]{\title{#1}\def\@ltitle{#1}}
\def\@ltitleshort{\@ltitle}
\newcommand{\titolcurt}[1]{\def\@ltitleshort{#1}}


\newcommand{\abstractcat}[1]{\def\l@abstractcat{#1}}
\newcommand{\rebut}[1]{\def\@rebut{#1}}
\newcommand{\acceptat}[1]{\def\@acceptat{#1}}
%Modifiquem el títol i fem la portada
\def\s@lt{\vskip 2mm plus 2mm}
\def\s@ltet{\vskip 1mm plus 1mm}

\def\@maketitle{%
%\setlength{\unitlength}{1mm}%
\thispagestyle{empty}

\begin{center} \LARGE
Universitat Polit\`ecnica de Catalunya\\[0.2em]
Facultat de Matem\`atiques i Estad\'{i}stica
\end{center}
\vspace{1cm}
\begin{center} \LARGE
Master in Advanced Mathematics and Mathematical Engineering\\[0.2em]
\huge Master's thesis
\end{center}
\vspace{2cm}
\begin{center}
{\color{blau}\Huge\bf\@title}\\[1em]
\LARGE \textbf{\l@authorStudent}
\end{center}
\vspace{4cm}
\begin{center}
\Large Supervised by \l@supervisors \\[0.5em]
\l@monthYear
\end{center}

\newpage\null
\thispagestyle{empty}

\newpage \null
\thispagestyle{empty}
\@graiments

\newpage \null \thispagestyle{empty} \newpage
\clearpage
\setcounter{page}{1}

\noindent
{\color{gris}\Large\bf Abstract}\\[1em]
\l@abstracteng\\[2cm]
{\color{gris}\large\bf Keywords}\\[1em]
\@paraulesclau

\newpage

}

% Guardem els autors a les variables
%
%\newcounter{num@utors}\setcounter{num@utors}{0}
%\newcommand{\autor}[1]{\author{#1}\addtocounter{num@utors}{1}\expandafter\def\csname @lautor\thenum@utors\endcsname{#1}}
%\newcommand{\univautor}[1]{\expandafter\def\csname l@univautor\thenum@utors\endcsname{#1}}
%\newcommand{\emailautor}[1]{\expandafter\def\csname l@emailautor\thenum@utors\endcsname{#1}}
%
\newcounter{num@utors}\setcounter{num@utors}{0}
\newcommand{\autor}[1]{\author{#1}\addtocounter{num@utors}{1}\expandafter\def\csname @lautor\thenum@utors\endcsname{#1}\expandafter\def\csname @pautor\thenum@utors\endcsname{#1}}
\newcommand{\corrautor}[1]{\author{#1}\addtocounter{num@utors}{1}\expandafter\def\csname @lautor\thenum@utors\endcsname{#1}\expandafter\def\csname @pautor\thenum@utors\endcsname{${}^\ast$#1}}
\newcommand{\univautor}[1]{\expandafter\def\csname l@univautor\thenum@utors\endcsname{#1}}
\newcommand{\emailautor}[1]{\expandafter\def\csname l@emailautor\thenum@utors\endcsname{#1}}
% Fem la caixa que apareixerà a les capçaleres
\newcounter{i@autor}
\newcommand{\autorsbox}{\setcounter{i@autor}{1}%
\loop
\bf\it\csname @lautor\thei@autor\endcsname % {\color{gris}\csname l@univautor\thei@autor\endcsname} \\
\ifnum \c@i@autor < \c@num@utors , \addtocounter{i@autor}{1} 
\repeat
}
% \newcommand{\autorsbox}{\setcounter{i@autor}{1}%
% \begin{minipage}[b]{15cm}
% \loop
% \bf\it\csname @lautor\thei@autor\endcsname,  {\color{gris}\csname l@univautor\thei@autor\endcsname} \\
% \ifnum \c@i@autor < \c@num@utors \addtocounter{i@autor}{1} 
% \repeat
% \end{minipage}
% }
%Modifiquem les capçaleres i els peus
\RequirePackage{fancyhdr}
%\newcommand{\authorsbox}{\begin{tabular}{ll}Albert Ruiz & Universitat Autònoma de Barcelona\\ Antonio Viruel & UMA\end{tabular}}
\newcommand{\titlebox}{\makebox[170mm][r]{\@title}}
\renewcommand{\headrulewidth}{0pt}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[RO]{\setlength{\unitlength}{1mm}%
\begin{picture}(0,0)
\put(-170,10){\autorsbox}
%\put(-190,0){\includegraphics[width=19cm]{../00estil/linia_blava_logo.png}}
\end{picture}}
\fancyhead[RE]{\setlength{\unitlength}{1mm}%
\begin{picture}(0,0)
\put(-170,10){\makebox[170mm][r]{\@ltitleshort}}
%\put(-190,7.5){\includegraphics[width=21cm,height=0.5mm]{../00estil/linia_blava.png}}
\end{picture}}
\fancyfoot[RO]{\setlength{\unitlength}{1mm}%
\begin{picture}(0,0)
%\put(-190,5){\includegraphics[width=21cm]{../00estil/linia_blava.png}}
\put(-20,-3){\makebox[20mm][r]{\it\thepage}}
%\put(-125,-3){\makebox[90mm][r]{\color{gris} Reports@SCM \textbf{\l@numero}\ (\l@any),\ \l@pagin--\l@pagfin;\ DOI:\l@doi.}}
\end{picture}}
\fancyfoot[RE]{\setlength{\unitlength}{1mm}%
\begin{picture}(0,0)
%\put(-170,5){\includegraphics[width=21cm]{../00estil/linia_blava.png}}
%\put(-130,-5){\includegraphics[height=6mm]{../00estil/SCMAmbText.jpg}}
%\put(-113,-5){\includegraphics[height=6mm]{../00estil/iec.jpg}}
%\put(-95,-3){\url{http://reportsascm.iec.cat}}
\put(-170,-3){\makebox[20mm][l]{\it\thepage}}
\end{picture}}
%Seccions i altres
\RequirePackage{titlesec}
\titleformat{\section}
          {\color{gris}\normalfont\huge\bfseries} 
          {\thesection.}
          {1ex \@plus.2ex}
          {}
\titleformat{\subsection}
          {\normalfont\Large\bfseries} 
          {\Large\thesubsection}
          {1ex \@plus.2ex}
          {}
\titleformat{\subsubsection}
          {\normalfont\large\bfseries} 
          {\huge\thesubsection}
          {1ex \@plus.2ex}
          {}
\titleformat{\paragraf}
          {\normalfont\normalsize\bfseries} 
          {\huge\thesubsection}
          {1ex \@plus.2ex}
          {}
\titleformat{\subparagraf}
          {\normalfont\normalsize\bfseries} 
          {\huge\thesubsection}
          {1ex \@plus.2ex}
          {}

% \renewcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
%                                     {3.25ex \@plus1ex \@minus.2ex}%
%                                     {-1em}%
%                                     {\normalfont\normalsize\bfseries}}
% \renewcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}%
%                                        {3.25ex \@plus1ex \@minus .2ex}%
%                                        {-1em}%                                     
%                                        {\normalfont\normalsize\bfseries}}
% No te indentacio i te separacio entre paragrafs
\setlength{\parskip}{1mm}
%\setlength{\parindent}{0pt}
% Modifiquem les fonts
\renewcommand*{\familydefault}{\sfdefault}
%\renewcommand{\normalsize}{\fontsize{12pt}{14pt}\selectfont}

% Modifiquem la bibliografia
\usepackage{multicol}
\usepackage{etoolbox}
%\patchcmd{\thebibliography}{\section*{\refname}}
%   {\begin{multicols}{2}[\section*{\refname}]}{}{}
%\patchcmd{\endthebibliography}{\endlist}{\endlist\end{multicols}}{}{}
%\renewcommand\@biblabel[1]{#1.}
%\makeatother
