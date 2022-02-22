<!--
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
-->

\page cs_dg_writing_theory Writing the theory documentation

[TOC]

General rules
=============

Theses general rules should be seen as basic golden rules helping the whole
documentation to be consistent. They are strongly recommended:

- Respect a plan where you first present a general overview of the theory
  (what is it about, what is the main goal), then you present the equations
  in general, and finally the specific choices you have made.
- Use the macros described in the [macros section](@ref pg_sec_theory_macros)`
  *i.e* `usepackage{csmacros}`).
- Use the notations defined in the nomenclature of the theory guide
  as much as possible.
- Focus on your specificities and cite the generalities (external to EDF!), which
  you should add to the `biblio.bib` file located in the `doc/style` directory.
- Write in English (UK for the theory manual, US for the rest,
  for consistence with existing documentation).
- Use the existing style of code_saturne, that is to say use the `csdoc.csl`
  class (for long documents as a report) `csshortdoc.cls` class
  (for short documents as an article).
- Respect \f$ \mbox{\LaTeX} \f$ philosophy, as it is designed to make sensible spacing
  decisions by itself, do not use explicit horizontal or vertical spacing
  commands, except in a few accepted (mostly mathematical) situations.
- keep your own macros to an absolute minimum.

Macros and typography{#pg_sec_theory_macros}
=====================

This section does not pretend to describe how to write a \f$ \mbox{\LaTeX} \f$ document,
but is to present the macros defined in `csmacro.sty` and give some typographic
pieces of advice.

Macros
------

The `\CS` macro in the `csdoc.sty` package is used to allow a short
syntax and typeset the code_saturne name in a proper and consistant manner.

The available macros for mathematical symbols are available through
the `csmacros.sty` package.

<table>
<caption id="latex_name_macro_op">Mathematical operators defined in csmacros</caption>
<tr><th> \f$ \mbox{\LaTeX} \f$ code  <th> preview <th> comment
<tr><td> `$\divs$`           <td> \f$ \divs    \f$ <td>
<tr><td> `$\divv$`           <td> \f$ \divv    \f$ <td>
<tr><td> `$\divt$`           <td> \f$ \divt    \f$ <td>
<tr><td> `$\grad$`           <td> \f$ \grad    \f$ <td>
<tr><td> `$\ggrad$`          <td> \f$ \ggrad   \f$ <td>
<tr><td> `$\gradv$`          <td> \f$ \gradv   \f$ <td>
<tr><td> `$\gradt$`          <td> \f$ \gradt   \f$ <td>
<tr><td> `$\gradtt$`         <td> \f$ \gradtt  \f$ <td>
<tr><td> `$\mat{M}$`         <td> \f$ \mat{M}  \f$ <td>
<tr><td> `$\matt{M}$`        <td> \f$ \matt{M} \f$ <td>
<tr><td> `$\rot$`            <td> \f$ \rot     \f$ <td>
<tr><td> `$\vect{V}$`        <td> \f$ \vect{V} \f$ <td>
<tr><td> `$\tens{T}$`        <td> \f$ \tens{T} \f$ <td>
<tr><td> `$\transpose{M}$`   <td> \f$ \transpose{M} \f$ <td>
<tr><td> `$\symmetric{M}$`   <td> \f$ \symmetric{M} \f$ <td>
<tr><td> `$\trace$`          <td> \f$ \trace   \f$ <td>
<tr><td> `$\deviator{M}$`    <td> \f$ \deviator{M}  \f$ <td>
<tr><td> `$\norm{M}$`        <td> \f$ \norm{M} \f$ <td>
<tr><td> `$\rans{M}$`        <td> \f$ \rans{M} \f$ <td>
<tr><td> `$\fluct{M}$`       <td> \f$ \fluct{M} \f$ <td>
<tr><td> `$\fluctt{M}$`      <td> \f$ \fluctt{M} \f$ <td>
<tr><td> `$\favre{M}$`       <td> \f$ \favre{M} \f$ <td>
<tr><td> `$\ints{M}{N}$`     <td> \f$ \int{M}{N} \f$ <td>
<tr><td> `$\intss{M}{N}$`    <td> \f$ \intss{M}{N} \f$ <td>
<tr><td> `$\intt{M}{N}$`     <td> \f$ \intt{M}{N} \f$ <td>
<tr><td> `\degresC`          <td> \f$\mbox{\degresC}\f$ <td>
<tr><td> `$\Max$`            <td> \f$ \Max     \f$ <td>
<tr><td> `$\Min$`            <td> \f$ \Min     \f$ <td>
<tr><td> `$\dd$`             <td> \f$ \dd      \f$ <td> total derivative
</table>

Many macros are dedicated to discretized quantity notations used throughout
code_saturne. The following table lists the main ones, but may not be complete,
so checking the actual contents of `csmacros.sty` is always recommened.

<table>
<caption id="latex_name_macro_q">Discretized quanties defined in csmacros</caption>
<tr><th> \f$ \mbox{\LaTeX} \f$ code          <th> preview <th> comment
<tr><td> `$\Facei{\celli}$` <td> \f$ \Facei{\celli} \f$ <td> set of internal faces
<tr><td> `$\Faceb{\cellj}$` <td> \f$ \Faceb{\cellj} \f$ <td> set of boundary faces
<tr><td> `$\Face{\celli}$`  <td> \f$ \Face{\celli}  \f$ <td> set of faces
<tr><td> `$\face$`          <td> \f$ \face   \f$ <td> face
<tr><td> `$\fij$`           <td> \f$ \fij    \f$ <td> internal face
<tr><td> `$\fib$`           <td> \f$ \fib    \f$ <td> boundary face
<tr><td> `$\iface$`         <td> \f$ \iface  \f$ <td> oriented face
<tr><td> `$\ij$`            <td> \f$ \ij     \f$ <td> oriented internal face
<tr><td> `$\ib$`            <td> \f$ \ib     \f$ <td> oriented boundary face
<tr><td> `$\celli$`         <td> \f$ \celli  \f$ <td> name of the current cell
<tr><td> `$\cellj$`         <td> \f$ \cellj  \f$ <td> name of the adjacent cell
<tr><td> `$\ipf$`           <td> \f$ \ipf    \f$ <td> orthogonal center index of the current cell
<tr><td> `$\jpf$`           <td> \f$ \jpf    \f$ <td> orthogonal center index of the adjacent cell
<tr><td> `$\centi$`         <td> \f$ \centi  \f$ <td> center of the current cell
<tr><td> `$\centj$`         <td> \f$ \centj  \f$ <td> center of the adjacent cell
<tr><td> `$\centip$`        <td> \f$ \centip \f$ <td> orthogonal center of the current cell
<tr><td> `$\centjp$`        <td> \f$ \centjp \f$ <td> orthogonal center of the adjacent cell
<tr><td> `$\cento$`         <td> \f$ \cento  \f$ <td> intersection between the cell centers and the face
<tr><td> `$\centf$`         <td> \f$ \centf  \f$ <td> center of the face
</table>

Typography
----------

Here are some useful tricks:

- If you want to describe multiple topics, use the
  `\begin{itemize} \item \end{itemize}` environment.
- You can use blue and orange EDF colors with the blue `\textcolor{blueedf}{text}`,
  its darker version `\textcolor{bluededf}{text}`, or the orange
  `\textcolor{orangeedf}{text}` and its the dark version
  `\textcolor{orangededf}{text}`.
- Use label and references, and dissociate equations with sections and appendices
  and figures and tables using `\label{eq:label}`, `\label{sec:label}`,
  `\label{ap:label}`, `\label{fig:label}` and `\label{tab:label}` prefixes.
- Use the `\emph{}` mode for acronyms (*e.g.*, *EDF*).
- Use the `\emph{}` mode for Latin words (*e.g.*, *i.e.*, *a priori*, *etc.*).
- Use `\left(` instead of `(` and `\right)` instead of `)` in math mode.
- **DO NOT** put a space before the ":"  symbol. In English the rule is no
  space, never.
- **DO NOT** use `\newline` or `\\` except in a tabular environment or an array.
- Write "Equation" with a  first upper case letter. Use `\figurename~` and
  `\tablename~` to write \f$ \figurename \f$ and \f$ \tablename \f$.
- Use the enumerate environment:
  ```{.tex}
  \begin{enumerate}[ label=\roman{*}/, ref=(\roman{*})]
  \item $1^{st}$ item
  \item $2^{nd}$ item
  \end{enumerate}

  ```
  i/ 1<sup>st</sup> item <br>
  ii/ 2<sup>nd</sup> item <br>

- Use the remarks `\begin{remark} \end{remark}`
  and example `\begin{example} \end{example}` environments defined in `csdoc.cls`:
  ```{.tex}
  \begin{remark}
    A remark
  \end{remark}
  \begin{example}
    An example
  \end{example}}
  ```
  **Remark 1.1** _A remark <br> 
  **Example 1.1** _An example_
