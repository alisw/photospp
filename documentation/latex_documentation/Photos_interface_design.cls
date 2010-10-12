\ProvidesClass{my_report}
\LoadClass[11pt, titlepage]{article}
\usepackage[paperwidth=21cm,paperheight=29.7cm,includehead,headheight=2.0cm,pdftex,hmargin={3cm,2.5cm},vmargin={0cm,2cm},]{geometry} 
\setlength{\parindent}{0cm} 
\renewcommand{\baselinestretch}{1.2}
\parskip=1em

\renewcommand{\contentsname}{Table of Contents}

\renewcommand{\normalsize}{\fontsize{11pt}{11pt}\selectfont}
\renewcommand{\title}[1]{\textsf{\Huge{#1}}}
\renewcommand{\author}[1]{\textsf{\Huge{#1}}}
\newcommand{\subtitle}[1]{\textsf{\textbf{\Large{#1}}}}
\newcommand{\abstractheading}[1]{\textsf{\textbf{\LARGE{#1}}}}
\newcommand{\code}[1]{\texttt{\footnotesize{#1}}}

\renewcommand{\section}{\@startsection
{section}%                    % the name
{1}%                          % the level
{0mm}%                        % the indent
{10.6mm}%             % the beforeskip
{4.2mm}%           % the afterskip
{\LARGE\bfseries\sffamily}}  % the style

\renewcommand{\subsection}{\@startsection
{subsection}%                    % the name
{2}%                          % the level
{0mm}%                        % the indent
{6.4mm}%             % the beforeskip
{1.1mm}%           % the afterskip
{\Large\bfseries\sffamily}}  % the style

\renewcommand{\subsubsection}{\@startsection
{subsubsection}%                    % the name
{3}%                          % the level
{0mm}%                        % the indent
{4.2mm}%             % the beforeskip
{1.1mm}%           % the afterskip
{\normalsize\bfseries\sffamily}}  % the style

%% Chapter headings should be centered, uppercase, and at the top of the page.
\def\@makechapterhead#1{%
  { \parindent \z@ \raggedright \normalfont
    %\centering
    \ifnum \c@secnumdepth >\m@ne
        \huge\textsf{\@chapapp\space \thechapter:}
        % \par\nobreak
        %\vskip 20\p@
    \fi
    \interlinepenalty\@M
    \huge \bfseries \textsf{#1}\par\nobreak
    \rule{5cm}{0.5pt}
    \vskip 20\p@
  } }
  
\def\@makeschapterhead#1{%
  %\vspace*{50\p@}%
  { \parindent \z@ \raggedright
    %\centering
    \normalfont
    \interlinepenalty\@M
    \huge \bfseries  \textsf{#1}\par\nobreak
    \rule{5cm}{0.5pt}
    \vskip 20\p@

  }}
  
 \renewenvironment{abstract}{%
      \chapter*{\abstractname}%
      \addcontentsline{toc}{chapter}{\abstractname}
 }
     
     
\makeatletter
\renewcommand{\l@section}{\bfseries\@dottedtocline{1}{0em}{2.3em}}
\renewcommand{\l@subsection}{\normalfont\@dottedtocline{2}{2em}{2.3em}}
\renewcommand{\l@subsubsection}{\normalfont\@dottedtocline{3}{2em}{2.3em}}
%\renewcommand{\l@subsubsection}{\normalfont\@dottedtocline{4}{2em}{2.3em}}
\makeatother

\def\maketitle{\begin{titlepage}
\thispagestyle{empty}
\let\footnotesize\small \let\footnoterule\relax \setcounter{page}{0}
\null
\vfil
\hskip 11.2 cm IFJPAN-IV-2010-6

\hskip 12 cm August 2010

\begin{center}
\title{PHOTOS Interface in C++\\[0.5cm]{} 
\textbf{\small Draft of Technical and Physics Documentation\\[0.8cm]\rule{4cm}{1pt}\\[0.5cm]}}

\author{\LARGE N. Davidson$^{a,b}$, T. Przedzinski$^c$, Z. Was$^a$}\\[0.5cm]
\textsf{\Large{\em $^a$ Institute of Nuclear Physics, 
Polish Academy of Sciences,\\
         ul. Radzikowskiego 152, 31-342 Cracow, Poland}\\
{\em $^b$ University of Melbourne, Department of Physics \\ Australia.  }\\
{\em $^c $ The Faculty of Physics, Astronomy and Applied Computer
Science,\\
Jagellonian University, Reymonta 4, 30-059 Cracow, Poland.}
 }\\[0.5cm]

%\includegraphics[height=3cm]{logo.eps}\\[1cm]
{\bf \Large Abstract} \vskip 2 mm 

\end{center}

{ The first version of {\tt PHOTOS} Monte Carlo for
bremsstrahlung in the decay of particles and resonances
 with  interface to the HepMC event record written in C++ is now available.
The current functionality is already as good as that of its FORTRAN 
predecessor. 
Better steering options are introduced. An algorithm for event record interface
is adopted 
to the forthcoming installation of process dependent weights 
necessary to implement complete first order matrix element effects.
A multitude of tests is collected. Brief presentation 
of physics assumptions is also given. In particular, it is briefly presented
how the program exploits properties of first and second order matrix elements 
and why the phase space parametrization used is exact. The main purpose of the 
paper is to present technical aspects of the {\tt PHOTOS} Monte Carlo use. }

\centerline{\normalsize \today}

\vfil

\footnoterule
\noindent
{\footnotesize
\begin{itemize}
\item[${\dag}$]

This work is partially supported by  EU Marie Curie Research Training Network 
grant under the contract No. MRTN-CT-2006-0355505 and by  
Polish Government grant N202 06434 (2008-2011). 

\end{itemize}
}
\null
\end{titlepage}}

\pagestyle{plain}

