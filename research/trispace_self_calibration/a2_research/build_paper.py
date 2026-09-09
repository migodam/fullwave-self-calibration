"""Render the controlled Markdown manuscript subset to editable LaTeX/PDF.

Requires tectonic. Scientific equations are copied verbatim; prose is escaped.
This builder does not alter research records or execute experiments.
"""
from pathlib import Path
import re
import subprocess
import json
import hashlib

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[2]
OUT=PROJECT/"output/pdf"
SOURCE=HERE/"PAPER_DRAFT_A2.md"

def escape(s):
    s=s.replace("\u2011","-").replace("\u2013","--").replace("\u2014","---")
    chars={"&":r"\&","%":r"\%","#":r"\#","_":r"\_","{":r"\{","}":r"\}","~":r"\textasciitilde{}","^":r"\textasciicircum{}","\\":r"\textbackslash{}"}
    return "".join(chars.get(c,c) for c in s)

def inline(s):
    tokens=re.split(r"(\$[^$]+\$|\[[^\]]+\]\(https?://[^)]+\)|`[^`]+`|\*\*[^*]+\*\*)",s)
    out=[]
    for t in tokens:
        if t.startswith("$") and t.endswith("$"): out.append(t)
        elif re.fullmatch(r"\[[^\]]+\]\(https?://[^)]+\)",t):
            m=re.fullmatch(r"\[([^\]]+)\]\(([^)]+)\)",t)
            out.append(r"\href{"+m[2]+"}{"+escape(m[1])+"}")
        elif t.startswith("`") and t.endswith("`"):out.append(r"\texttt{"+escape(t[1:-1])+"}")
        elif t.startswith("**") and t.endswith("**"):out.append(r"\textbf{"+escape(t[2:-2])+"}")
        else:out.append(escape(t))
    return "".join(out)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    source=SOURCE.read_text();lines=source.splitlines();title=lines[0][2:]
    assert "PENDING" not in source
    out=[r"\documentclass[11pt,a4paper]{article}",r"\usepackage[margin=22mm]{geometry}",r"\usepackage{amsmath,amssymb,mathtools}",r"\usepackage{graphicx,booktabs,tabularx,array}",r"\usepackage{fontspec}",r"\setmainfont{TeX Gyre Termes}",r"\usepackage{xurl,hyperref}",r"\hypersetup{colorlinks=true,urlcolor=blue,linkcolor=black,pdftitle={"+escape(title)+r"}}",r"\usepackage{caption,enumitem,needspace,fancyhdr}",r"\captionsetup{font=small,labelfont=bf}",r"\setlist{nosep,leftmargin=*}",r"\setlength{\parindent}{0pt}",r"\setlength{\parskip}{5pt}",r"\setlength{\emergencystretch}{2em}",r"\setlength{\headheight}{15pt}",r"\pagestyle{fancy}\fancyhf{}\fancyhead[L]{\small TriSpace SOM / full-wave self-calibration}\fancyhead[R]{\small Research draft}\fancyfoot[C]{\thepage}",r"\renewcommand{\headrulewidth}{0.3pt}",r"\newsavebox{\eqbox}",r"\newcommand{\fitdisplay}[1]{\sbox{\eqbox}{$\displaystyle #1$}\ifdim\wd\eqbox>\linewidth\resizebox{\linewidth}{!}{\usebox{\eqbox}}\else\usebox{\eqbox}\fi}",r"\begin{document}",r"\begin{center}{\LARGE\bfseries "+escape(title)+r"\par}\vspace{8pt}{\small 6 September 2026 | Audited research draft}\end{center}"]
    out=[x.replace("TeX Gyre Termes","Times New Roman") for x in out]
    i=1;tables=0;figures=0;equations=0
    while i<len(lines):
        line=lines[i].strip()
        if not line:i+=1;continue
        if line=="$$":
            end=lines.index("$$",i+1);eq="\n".join(lines[i+1:end]);equations+=1
            out.append("\\begin{equation}\n\\fitdisplay{"+eq+"}\n\\end{equation}");i=end+1;continue
        if line.startswith("## ") or line.startswith("### "):
            level="subsection" if line.startswith("### ") else "section"
            txt=line.lstrip("# ")
            if txt.startswith("Appendix A"):
                out.append(r"\clearpage")
            out.append(r"\needspace{5\baselineskip}"+"\\"+level+"*{"+inline(txt)+"}");i+=1;continue
        image=re.fullmatch(r"!\[(.*)\]\((.*)\)",line)
        if image:
            path=(HERE/image[2]).resolve();assert path.exists(),path
            out += [r"\begin{figure}[htbp]\centering",r"\includegraphics[width=\linewidth,height=.65\textheight,keepaspectratio]{\detokenize{"+str(path)+"}}",r"\caption{"+inline(image[1])+r"}\end{figure}"];figures+=1;i+=1;continue
        if line.startswith("|"):
            rows=[]
            while i<len(lines) and lines[i].strip().startswith("|"):
                row=[c.strip() for c in lines[i].strip().strip("|").split("|")]
                if not all(re.fullmatch(r":?-+:?",c) for c in row):rows.append(row)
                i+=1
            n=len(rows[0]);assert all(len(r)==n for r in rows)
            tables+=1
            # Six-column budget table gives method names additional width.
            if n==6:
                widths=[.55,1.8,.8,.95,1.,.9]
                fmt="".join(r">{\hsize="+str(w)+r"\hsize\linewidth=\hsize\raggedright\arraybackslash}X" for w in widths)
            else:fmt=r">{\raggedright\arraybackslash}X"*n
            out += [r"\begin{table}[htbp]\centering\small",r"\setlength{\tabcolsep}{4pt}\renewcommand{\arraystretch}{1.2}",r"\begin{tabularx}{\linewidth}{"+fmt+r"}\toprule"]
            for j,row in enumerate(rows):
                vals=[inline(c) for c in row]
                if j==0:vals=[r"\textbf{"+v+"}" for v in vals]
                out.append(" & ".join(vals)+r" \\")
                if j==0:out.append(r"\midrule")
            out += [r"\bottomrule\end{tabularx}\end{table}"]
            continue
        if re.match(r"\d+\. ",line):
            out.append(r"\begin{enumerate}")
            while i<len(lines) and re.match(r"\d+\. ",lines[i].strip()):
                out.append(r"\item "+inline(re.sub(r"^\d+\. ","",lines[i].strip())));i+=1
            out.append(r"\end{enumerate}");continue
        out.append(inline(line)+"\n");i+=1
    out.append(r"\end{document}")
    tex=OUT/"trispace_a2_research_draft.tex"
    tex.write_text("\n".join(out)+"\n")
    subprocess.run(["/opt/homebrew/bin/tectonic","--keep-logs","--outdir",str(OUT),str(tex)],check=True)
    meta=dict(source=str(SOURCE),source_sha256=hashlib.sha256(source.encode()).hexdigest(),equations=equations,tables=tables,figures=figures,pdf=str(tex.with_suffix(".pdf")))
    (OUT/"build_manifest.json").write_text(json.dumps(meta,indent=2)+"\n")
    print(json.dumps(meta,indent=2))

if __name__=="__main__":main()
