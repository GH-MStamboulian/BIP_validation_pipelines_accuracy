\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{protocol}[Guardant class for study protocols]
 
% Extend the LaTeX articl class
\LoadClass[a4paper,12pt]{article}

% Set the whitespace for margins, headers, and footers
\RequirePackage{geometry}
\geometry{letterpaper,
          left=0.5in,
          right=0.5in,
          top=2in,
          headheight=1.25in,
          headsep=0.25in,
          bmargin=0.75in}

% Make all of the hyperlinks black
\RequirePackage[colorlinks=black,
                urlcolor=black,
                linkcolor=black,
                citecolor=black,
                unicode]{hyperref}
                
% Recreated the header from Guardants Microsoft Word documents
\RequirePackage{lastpage}
\def\headertitle#1{\def\@headertitle{#1}}
\def\documentid#1{\def\@documentid{#1}}
\def\revisionnumber#1{\def\@revisionnumber{#1}}
\RequirePackage{fancyhdr}
\pagestyle{fancy}
\renewcommand{\headrulewidth}{0pt}
\fancyhead[CO,LO,RO]{} %% clear out all headers
\fancyhead[C]{%
          \begin{tabular}{|C{2.2in}|C{3.1in}|C{1.7in}|C{0.1in}}
          \cline{1-3}
          \includegraphics{logo.png} &
          \multicolumn{2}{C{4.8in}|}{\textbf{\large \@headertitle}} &
          \vspace{0.6in} \\
          \cline{1-3}
          \footnotesize \textbf{Document ID:\newline \@documentid} &
          \footnotesize \textbf{Revision Number:\newline \@revisionnumber} &
          \footnotesize \textbf{Page:\newline \thepage\ of \pageref{LastPage}} &
          \vspace{0.15in} \\[0.1in]
          \cline{1-3}
          \end{tabular}%
}
\fancyfoot[CO,LO,RO]{} %% clear out all footers
\fancyfoot[C]{%
    \textbf{Confidential Property of Guardant Health, Inc.}
}

% Command to reference an appendix with number as well as its title
\newcommand{\appendixref}[1]{Appendix \autoref{#1}. \nameref{#1}}
\newcommand{\shortappendixref}[1]{\textbf{Appendix \autoref{#1}}}

 % Add a blank line before each paragraph
\setlength{\parskip}{1em}
\setlength{\parindent}{0em}

% Begin enumerated lists according to the section number
\RequirePackage{enumitem}
\setenumerate[1]{label=\thesection.\arabic*.}
\setenumerate[2]{label*=\arabic*.}
\setenumerate[3]{label*=\arabic*.}

% Make a column type for tables that is fixed width and centered
\RequirePackage{array}
\newcolumntype{C}[1]{>{\centering\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}
\newcolumntype{L}[1]{>{\raggedright\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}

% Change the spacing before and after section headers
\RequirePackage{titlesec}
\titlespacing\section{0pt}{5pt plus 4pt minus 2pt}{-6pt plus 2pt minus 2pt}
% Change the format of the section titles
\titlelabel{\thetitle.\quad}
\titleformat{\section}
  {}{\textbf{\thesection.}}{3pt}{\textbf}
% and also for subsection titles
\titleformat{\subsection}
  {}{\textbf{}}{3pt}{\textbf}
\titleformat{\subsubsection}
  {}{\textbf{}}{3pt}{\textbf}
\titlespacing\subsection{0pt}{5pt plus 4pt minus 2pt}{-3pt plus 2pt minus 2pt}
\titlespacing\subsubsection{0pt}{5pt plus 4pt minus 2pt}{0pt plus 2pt minus 2pt}

% Better cross referencing
\RequirePackage{amsmath}
\RequirePackage[capitalise]{cleveref}
\crefdefaultlabelformat{#2\textbf{#1}#3}
\creflabelformat{equation}{#2\textup{\textbf{#1}}#3}
\crefname{section}{\textbf{Section}}{\textbf{Sections}}
\Crefname{section}{\textbf{Section}}{\textbf{Sections}}
\crefname{table}{\textbf{Table}}{\textbf{Tables}}
\Crefname{table}{\textbf{Table}}{\textbf{Tables}}
\crefname{appsec}{\textbf{Appendix}}{\textbf{Appendices}}
\Crefname{appsec}{\textbf{Appendix}}{\textbf{Appendices}}
\crefname{equation}{\textbf{Equation}}{\textbf{Equations}}
\Crefname{equation}{\textbf{Equation}}{\textbf{Equations}}
\crefname{enumi}{\textbf{Section}}{\textbf{Sections}}
\Crefname{enumi}{\textbf{Section}}{\textbf{Sections}}
\crefname{figure}{\textbf{Figure}}{\textbf{Figures}}
\Crefname{figure}{\textbf{Figure}}{\textbf{Figures}}

% Allow for different fonts
\RequirePackage{fontspec}
