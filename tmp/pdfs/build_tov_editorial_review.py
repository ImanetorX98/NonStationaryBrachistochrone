#!/usr/bin/env python3

from __future__ import annotations

import os
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    KeepTogether,
    PageBreak,
    PageTemplate,
    Paragraph,
    Preformatted,
    Spacer,
    Table,
    TableStyle,
)


OUTPUT_PATH = os.path.abspath(
    "TOVInternalMetric/output/pdf/TOVBrach_editorial_review_proposed_changes.pdf"
)
STIX_REGULAR = "/System/Library/Fonts/Supplemental/STIXTwoText.ttf"
STIX_ITALIC = "/System/Library/Fonts/Supplemental/STIXTwoText-Italic.ttf"


def register_fonts() -> None:
    pdfmetrics.registerFont(TTFont("STIXText", STIX_REGULAR))
    pdfmetrics.registerFont(TTFont("STIXTextItalic", STIX_ITALIC))


def page_footer(canvas_obj, document) -> None:
    canvas_obj.saveState()
    width, _ = A4
    canvas_obj.setStrokeColor(colors.HexColor("#B8BDC6"))
    canvas_obj.setLineWidth(0.45)
    canvas_obj.line(18 * mm, 14 * mm, width - 18 * mm, 14 * mm)
    canvas_obj.setFont("Helvetica", 7.5)
    canvas_obj.setFillColor(colors.HexColor("#5B6472"))
    canvas_obj.drawString(18 * mm, 9.5 * mm, "Editorial review - TOV brachistochrone manuscript")
    canvas_obj.drawRightString(width - 18 * mm, 9.5 * mm, f"Page {document.page}")
    canvas_obj.restoreState()


def build_styles():
    styles = getSampleStyleSheet()
    styles.add(
        ParagraphStyle(
            name="ReviewTitle",
            parent=styles["Title"],
            fontName="Helvetica-Bold",
            fontSize=21,
            leading=25,
            textColor=colors.HexColor("#12233F"),
            alignment=TA_LEFT,
            spaceAfter=10,
        )
    )
    styles.add(
        ParagraphStyle(
            name="ReviewSubtitle",
            parent=styles["Normal"],
            fontName="STIXText",
            fontSize=11.2,
            leading=15,
            textColor=colors.HexColor("#465366"),
            spaceAfter=14,
        )
    )
    styles.add(
        ParagraphStyle(
            name="Section",
            parent=styles["Heading1"],
            fontName="Helvetica-Bold",
            fontSize=14,
            leading=17,
            textColor=colors.HexColor("#12233F"),
            spaceBefore=12,
            spaceAfter=7,
            keepWithNext=True,
        )
    )
    styles.add(
        ParagraphStyle(
            name="Subsection",
            parent=styles["Heading2"],
            fontName="Helvetica-Bold",
            fontSize=10.8,
            leading=14,
            textColor=colors.HexColor("#27476E"),
            spaceBefore=9,
            spaceAfter=4,
            keepWithNext=True,
        )
    )
    styles.add(
        ParagraphStyle(
            name="BodySTIX",
            parent=styles["BodyText"],
            fontName="STIXText",
            fontSize=9.4,
            leading=13.1,
            textColor=colors.HexColor("#20242B"),
            alignment=TA_LEFT,
            spaceAfter=5,
        )
    )
    styles.add(
        ParagraphStyle(
            name="SmallSTIX",
            parent=styles["BodyText"],
            fontName="STIXText",
            fontSize=8.2,
            leading=11.2,
            textColor=colors.HexColor("#39414D"),
            spaceAfter=3,
        )
    )
    styles.add(
        ParagraphStyle(
            name="Formula",
            parent=styles["BodyText"],
            fontName="STIXText",
            fontSize=11.0,
            leading=15,
            alignment=TA_CENTER,
            textColor=colors.HexColor("#111827"),
            leftIndent=8 * mm,
            rightIndent=8 * mm,
            spaceBefore=4,
            spaceAfter=6,
        )
    )
    styles.add(
        ParagraphStyle(
            name="ReviewCode",
            fontName="Courier",
            fontSize=7.25,
            leading=9.5,
            textColor=colors.HexColor("#172033"),
            leftIndent=4 * mm,
            rightIndent=4 * mm,
            spaceBefore=4,
            spaceAfter=7,
            backColor=colors.HexColor("#F4F6F8"),
            borderColor=colors.HexColor("#D7DCE2"),
            borderWidth=0.5,
            borderPadding=6,
        )
    )
    styles.add(
        ParagraphStyle(
            name="Checklist",
            parent=styles["BodyText"],
            fontName="STIXText",
            fontSize=9.2,
            leading=12.5,
            leftIndent=5 * mm,
            firstLineIndent=-4 * mm,
            bulletIndent=0,
            spaceAfter=3,
        )
    )
    return styles


def paragraph(text: str, styles, style_name: str = "BodySTIX") -> Paragraph:
    return Paragraph(text, styles[style_name])


def bullet(text: str, styles) -> Paragraph:
    return Paragraph(f"&#8226;&nbsp; {text}", styles["Checklist"])


def callout(title: str, body: str, styles, severity: str = "major") -> Table:
    palette = {
        "critical": ("#8F1D1D", "#FFF1F0", "#F2B8B5"),
        "major": ("#8A4B08", "#FFF7E6", "#F2CF91"),
        "minor": ("#245A78", "#EEF7FB", "#A9D3E5"),
        "ok": ("#27623A", "#EEF8F1", "#ACD6B7"),
    }
    title_color, background, border = palette[severity]
    contents = [
        Paragraph(f"<b><font color='{title_color}'>{escape(title)}</font></b>", styles["BodySTIX"]),
        Paragraph(body, styles["SmallSTIX"]),
    ]
    table = Table([[contents]], colWidths=[174 * mm])
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(background)),
                ("BOX", (0, 0), (-1, -1), 0.7, colors.HexColor(border)),
                ("LEFTPADDING", (0, 0), (-1, -1), 8),
                ("RIGHTPADDING", (0, 0), (-1, -1), 8),
                ("TOPPADDING", (0, 0), (-1, -1), 6),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
            ]
        )
    )
    return table


def code_block(text: str, styles) -> Preformatted:
    return Preformatted(text.strip("\n"), styles["ReviewCode"], maxLineLength=105)


def build_story(styles):
    story = []
    story.append(Spacer(1, 8 * mm))
    story.append(paragraph("Editorial review and proposed changes", styles, "ReviewTitle"))
    story.append(
        paragraph(
            "Coordinate-Time and Proper-Time Brachistochrones in Tolman-Oppenheimer-Volkoff Stellar Interiors",
            styles,
            "ReviewSubtitle",
        )
    )
    story.append(
        callout(
            "Overall assessment",
            "The manuscript is scientifically coherent and substantially clearer than the earlier versions. "
            "The first-order t-tau splitting formula and its constant-density limit are correct. Before submission, "
            "the manuscript should nevertheless receive two mathematical-formal corrections and a small set of "
            "typesetting and exposition revisions. The sign typo after Eq. (74), the justification of the moving "
            "turning-point limit, the overlap in Eq. (67), and the visible hyperlink borders are the highest-priority items.",
            styles,
            "major",
        )
    )
    story.append(Spacer(1, 5 * mm))
    story.append(paragraph("Scope of this document", styles, "Section"))
    story.append(
        paragraph(
            "This report is written as an implementation brief. It identifies the exact source locations in "
            "<font name='Courier'>TOVBrachistocrone_copia.tex</font>, supplies replacement text or LaTeX where useful, "
            "and separates required corrections from optional editorial improvements. Line numbers refer to the "
            "August 4, 2026 source and will shift after editing.",
            styles,
        )
    )
    story.append(paragraph("Priority summary", styles, "Section"))
    priority_data = [
        ["Priority", "Action", "Reason"],
        ["P0", "Correct the sign after Eq. (74)", "The explanatory identity currently contradicts W = Psi - 1."],
        ["P0", "Regularize before differentiating F(q)", "The moving lower limit has an integrable square-root singularity."],
        ["P0", "Repair Eq. (67) layout", "The equation visibly overlaps the adjacent column."],
        ["P0", "Hide hyperlink borders", "Red, green, and cyan PDF rectangles are visible throughout."],
        ["P1", "Qualify the universality claim", "The result assumes a regular Newtonian limit and nondegenerate root."],
        ["P1", "Clarify proper-time minimization", "Avoid confusion with proper-time maximization between fixed events."],
        ["P1", "Standardize compactness mu", "The code contains both mu = M/R and mu = 2M/R conventions."],
        ["P2", "Reorder and condense sections", "The main analytic result arrives late; Expected Observables is too short."],
    ]
    priority_table = Table(priority_data, colWidths=[15 * mm, 58 * mm, 101 * mm], repeatRows=1)
    priority_table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#183B63")),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("FONTNAME", (0, 1), (-1, -1), "STIXText"),
                ("FONTSIZE", (0, 0), (-1, -1), 7.6),
                ("LEADING", (0, 0), (-1, -1), 10),
                ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#C9CED6")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F7F8FA")]),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ]
        )
    )
    story.append(priority_table)

    story.append(PageBreak())
    story.append(paragraph("1. Required mathematical corrections", styles, "Section"))
    story.append(paragraph("1.1 Sign after Eq. (74)", styles, "Subsection"))
    story.append(
        callout(
            "Required correction - source line 1266",
            "The manuscript currently states Psi(q<sub>0</sub>) - Psi(x) = W(x) - W(q<sub>0</sub>). "
            "Because W = Psi - 1, the right-hand side must have the opposite order. Equations (73), (74), and "
            "(76)-(78) already use the correct sign implicitly, so the central result is unchanged.",
            styles,
            "critical",
        )
    )
    story.append(paragraph("Current text:", styles, "SmallSTIX"))
    story.append(code_block(r"$\Psi(q_0)-\Psi(x)=W(x)-W(q_0)$.", styles))
    story.append(paragraph("Replace with:", styles, "SmallSTIX"))
    story.append(code_block(r"$\Psi(q_0)-\Psi(x)=W(q_0)-W(x)$.", styles))
    story.append(
        paragraph(
            "For x &gt; q<sub>0</sub> and positive density, W decreases outward, so W(q<sub>0</sub>) - W(x) is "
            "non-negative. This is consistent with the positive constant-density integrand in Eq. (76).",
            styles,
        )
    )

    story.append(paragraph("1.2 Differentiation with a moving singular endpoint", styles, "Subsection"))
    story.append(
        callout(
            "Required strengthening - source lines 1195-1200",
            "The conclusion is correct, but integrability alone does not justify differentiating an integral whose "
            "lower limit is also a square-root singular point. Introduce the same endpoint-regularizing substitution "
            "used in the numerical section before taking the derivative with respect to q.",
            styles,
            "major",
        )
    )
    story.append(paragraph("Suggested replacement paragraph:", styles, "SmallSTIX"))
    story.append(
        code_block(
            r"""
To justify the expansion with a moving turning point, first introduce
$x=q+(1-q)\sin^2 u$, with $u\in[0,\pi/2]$.  In the transformed
integral the interval is fixed and the square-root endpoint singularity
is removed.  The regularized functional is differentiable with respect
to $q$, and therefore
\[
\mathcal F_i(q_0+\mu q_1^{(i)},\mu)
=\mathcal F^{(0)}(q_0)
+\mu q_1^{(i)}\partial_q\mathcal F^{(0)}(q_0)
+\mu\!\int_{q_0}^{1}\!\mathcal I(x;q_0)\Theta_i(x;q_0)\,dx
+O(\mu^2).
\]
Thus the dependence generated by the displacement of the turning point
is contained in $\partial_q\mathcal F^{(0)}(q_0)$, while the explicit
first-order correction is evaluated at the fixed leading-order root $q_0$.
""",
            styles,
        )
    )
    story.append(
        paragraph(
            "Also state explicitly that the selected branch satisfies "
            "<font name='STIXTextItalic'>partial</font><sub>q</sub> F<super>(0)</super>(q<sub>0</sub>) != 0. "
            "Without this nondegeneracy condition, perturbative inversion is not defined.",
            styles,
        )
    )

    story.append(paragraph("1.3 Transparent derivation of delta_t - delta_tau", styles, "Subsection"))
    story.append(
        paragraph(
            "Equation (63) is correct. A two-line expansion would make the cancellation of the unknown second-order "
            "lapse contribution completely explicit. This is useful because Eq. (63) is the conceptual reason why the "
            "splitting can be computed without the full first post-Newtonian pressure correction.",
            styles,
        )
    )
    story.append(paragraph("Suggested insertion after Eq. (63):", styles, "SmallSTIX"))
    story.append(
        code_block(
            r"""
More explicitly, write
$A_0-A=\mu W\,[1+\mu\chi(x)+O(\mu^2)]$.  Then
\[
\delta_t=\frac{\Psi}{2}-\frac{\chi}{2},\qquad
\delta_\tau=-\frac{\Psi}{2}-\frac{\chi}{2},
\]
so that the common, pressure-dependent correction $\chi(x)$ cancels and
$\delta_t-\delta_\tau=\Psi$.
""",
            styles,
        )
    )

    story.append(PageBreak())
    story.append(paragraph("2. Verification of the central result", styles, "Section"))
    story.append(
        callout(
            "No change to Eqs. (75) and (78)",
            "The functional form, overall sign, and constant-density coefficient are correct. The issue is confined "
            "to the explanatory sign typo after Eq. (74) and to the rigor of the endpoint differentiation argument.",
            styles,
            "ok",
        )
    )
    story.append(paragraph("The verified chain is:", styles, "BodySTIX"))
    verification_items = [
        "A = 1 - mu Psi + O(mu<super>2</super>) and A<sub>0</sub> - A = mu W + O(mu<super>2</super>), with W = Psi - 1.",
        "n<sub>t</sub>/n<sub>tau</sub> = 1/A implies delta<sub>t</sub> - delta<sub>tau</sub> = Psi.",
        "The common square-root B contribution cancels in q<sub>1</sub><super>(t)</super> - q<sub>1</sub><super>(tau)</super>.",
        "Theta<sub>t</sub> - Theta<sub>tau</sub> = P<sub>0</sub><super>2</super>(x)[Psi(q<sub>0</sub>) - Psi(x)]/[P<sub>0</sub><super>2</super>(x) - P<sub>0</sub><super>2</super>(q<sub>0</sub>)].",
        "Therefore q<super>(t)</super> - q<super>(tau)</super> = -mu G[q<sub>0</sub>]/F<super>(0)</super>'(q<sub>0</sub>) + O(mu<super>2</super>).",
    ]
    for verification_item in verification_items:
        story.append(bullet(verification_item, styles))
    story.append(
        paragraph(
            "For constant density, W(x) = (1 - x<super>2</super>)/2, "
            "F<super>(0)</super>(q) = pi(1 - q)/2, and",
            styles,
        )
    )
    story.append(
        paragraph(
            "G[q<sub>0</sub>] = pi q<sub>0</sub>(1 - q<sub>0</sub><super>2</super>)<super>2</super>/8,",
            styles,
            "Formula",
        )
    )
    story.append(
        paragraph(
            "which gives q<super>(t)</super> - q<super>(tau)</super> = "
            "mu q<sub>0</sub>(1 - q<sub>0</sub><super>2</super>)<super>2</super>/4 + O(mu<super>2</super>). "
            "At Delta = pi/2, q<sub>0</sub> = 1/2 and the coefficient is 9/128.",
            styles,
        )
    )

    story.append(paragraph("3. Conditions behind the word 'universal'", styles, "Section"))
    story.append(
        paragraph(
            "The claim is best described as EOS-agnostic at fixed Newtonian density profile, rather than literally "
            "valid for any imaginable equation of state. The derivation assumes:",
            styles,
        )
    )
    for condition in [
        "a static, spherical configuration admitting a regular Newtonian limit;",
        "p/(rho c<super>2</super>) = O(mu), so pressure first affects the lapse beyond the leading coefficient a<sub>1</sub>;",
        "a smooth, non-negative normalized Newtonian density profile;",
        "a unique regular turning-point branch q<sub>0</sub>(Delta);",
        "a nonzero derivative F<super>(0)</super>'(q<sub>0</sub>).",
    ]:
        story.append(bullet(condition, styles))
    story.append(paragraph("Suggested abstract wording:", styles, "SmallSTIX"))
    story.append(
        code_block(
            r"""
A weak-field expansion in $\mu=r_s/R$ shows that, for any regular
barotropic stellar sequence admitting a smooth Newtonian limit, the
leading difference $q^{(t)}-q^{(\tau)}$ is determined by the normalized
Newtonian density profile through two one-dimensional quadratures and
one derivative, without integrating the full TOV system.
""",
            styles,
        )
    )
    story.append(paragraph("Suggested opening of Sec. X:", styles, "SmallSTIX"))
    story.append(
        code_block(
            r"""
We assume a regular weak-field stellar sequence for which
$p/(\rho c^2)=O(\mu)$ and the normalized Newtonian density profile
$\bar\rho(x)$ has a smooth $\mu\to0$ limit.  We further restrict to a
nondegenerate turning-point branch satisfying
$\partial_q\mathcal F^{(0)}(q_0)\ne0$.
""",
            styles,
        )
    )
    story.append(
        paragraph(
            "Define bar-rho explicitly as the Newtonian rest-mass density profile evaluated at r = Rx, up to an "
            "arbitrary normalization that cancels in bar-m. Replace the phrase 'expanded to O(mu<super>0</super>) in "
            "the pressure' with an explicit ordering statement: the pressure term is smaller than the mass-density "
            "term by O(mu) and therefore first contributes to A at O(mu<super>2</super>).",
            styles,
        )
    )

    story.append(PageBreak())
    story.append(paragraph("4. Physical exposition of t versus tau", styles, "Section"))
    story.append(
        paragraph(
            "Sections III-VI explain the distinction well: t is the Schwarzschild time normalized at infinity, whereas "
            "tau is the particle's worldline time. One likely referee objection should nevertheless be anticipated. "
            "Timelike geodesics maximize proper time between fixed spacetime events, while the paper minimizes tau. "
            "There is no contradiction because the endpoints are fixed only spatially and different tunnels produce "
            "different arrival events; moreover the tunnel supplies normal constraint forces, so the worldline is not "
            "an unconstrained spacetime geodesic.",
            styles,
        )
    )
    story.append(paragraph("Suggested insertion in Sec. V after Eq. (34):", styles, "SmallSTIX"))
    story.append(
        code_block(
            r"""
This variational problem should not be confused with the standard
proper-time extremization between two fixed spacetime events.  Here the
endpoints are fixed only on the stellar surface, while the arrival event
depends on the chosen tunnel; in addition, the tunnel exerts normal
constraint forces.  The usual geodesic proper-time maximization theorem
therefore does not apply to this comparison.
""",
            styles,
        )
    )
    story.append(
        paragraph(
            "A compact comparison table after Eq. (43) would further help the reader: minimized quantity, effective "
            "index, clock interpretation, and relation n<sub>t</sub>/n<sub>tau</sub> = 1/A.",
            styles,
        )
    )

    story.append(paragraph("5. Typesetting corrections", styles, "Section"))
    story.append(paragraph("5.1 Eq. (67) overlap", styles, "Subsection"))
    story.append(
        paragraph(
            "The current one-line equation exceeds the REVTeX column width and overlaps the text in the right column. "
            "Replace the equation environment at source lines 1164-1173 with:",
            styles,
        )
    )
    story.append(
        code_block(
            r"""
\begin{equation}
\begin{aligned}
\frac{1}{\sqrt{P_i^2(x)-P_i^2(q)}}
&=\frac{1}{\sqrt{D_0}}
\left[1-\mu\,
\frac{P_0^2(x)\delta_i(x)-P_0^2(q)\delta_i(q)}{D_0}
\right]
\\
&\quad +O(\mu^2).
\end{aligned}
\label{eq:denom_inv}
\end{equation}
""",
            styles,
        )
    )
    story.append(paragraph("5.2 Hyperlink borders", styles, "Subsection"))
    story.append(
        paragraph(
            "The PDF visibly prints red boxes around equation references, green boxes around citations, and cyan boxes "
            "around URLs and the email address. Replace the preamble line",
            styles,
        )
    )
    story.append(code_block(r"\usepackage{hyperref}", styles))
    story.append(paragraph("with either", styles, "SmallSTIX"))
    story.append(code_block(r"\usepackage[hidelinks]{hyperref}", styles))
    story.append(paragraph("or add", styles, "SmallSTIX"))
    story.append(code_block(r"\hypersetup{hidelinks}", styles))
    story.append(paragraph("5.3 Figure typography", styles, "Subsection"))
    story.append(
        paragraph(
            "The figures are legible, but their sans-serif typography differs from the REVTeX body. For a more uniform "
            "APS presentation, use a serif/STIX math font or configure Matplotlib to match the manuscript. This is an "
            "editorial improvement rather than a correctness issue.",
            styles,
        )
    )

    story.append(PageBreak())
    story.append(paragraph("6. Structure and captions", styles, "Section"))
    structural_actions = [
        "Move the universal weak-field section before the anisotropic application. The introduction presents the splitting formula as the main result, but the reader reaches it only after a long specialized detour.",
        "Merge the short Expected Observables section into the opening of the Discussion. It currently reads as a list rather than an independent scientific section.",
        "In the Fig. 1 caption, write that q^(t) > q^(tau) is observed for the displayed model. Avoid suggesting that pointwise n_t > n_tau alone proves the turning-radius ordering; the Discussion correctly notes that no general comparison theorem has been established.",
        "Revise the Fig. 7 caption. The Gamma -> infinity convergence validates the incompressible numerical limit, not by itself the weak-field coefficient in Eq. (78). The latter is tested by the small-mu slope in Fig. 6.",
        "If Eq. (75) is invoked in the anisotropic section, add a sentence explaining why anisotropic pressure does not enter the leading lapse coefficient, or describe the weak anisotropy dependence as a numerical observation rather than a direct prediction of the isotropic derivation.",
        "Consider limiting the visual emphasis of Fig. 6 to mu approximately below 0.2, or shade the weak-field region. At larger compactness the linear approximation visibly fails, which is expected but should not dominate the presentation of a first-order result.",
    ]
    for structural_action in structural_actions:
        story.append(bullet(structural_action, styles))

    story.append(paragraph("Suggested Fig. 7 final sentence:", styles, "SmallSTIX"))
    story.append(
        code_block(
            r"""
The convergence of the polytropic curves toward the incompressible
Schwarzschild solution checks the numerical constant-density limit; the
weak-field coefficient in Eq. (78) is tested independently by the
small-$\mu$ comparison in Fig. 6.
""",
            styles,
        )
    )

    story.append(paragraph("7. Notation and reproducibility", styles, "Section"))
    story.append(
        callout(
            "Compactness convention",
            "The manuscript consistently defines mu = r_s/R = 2GM/(Rc^2). The generator still contains legacy "
            "figures labelled mu = M/R as well as newer paper figures using mu = 2M/R. Internally retain "
            "M_over_R = M/R and reserve the name mu exclusively for 2*M_over_R. Otherwise a plot can be generated "
            "with an x-axis differing by a factor of two.",
            styles,
            "major",
        )
    )
    for reproducibility_item in [
        "The public tag v1.0-submission points to commit 5507e4b.",
        "The remote genera_grafici_tov.py at that commit is byte-identical to the local source inspected for this review.",
        "Manually open the newly minted Zenodo DOI once more before submission; the automated review environment could not independently resolve the fresh record.",
        "Where dimensional units are retained in the surrounding text, prefer ct/R and c tau/R over t/R and tau/R, even if a table caption states G = c = 1.",
        "Confirm the exact model names in the AI-use disclosure and align the wording with the current policy of the target journal.",
    ]:
        story.append(bullet(reproducibility_item, styles))

    story.append(paragraph("8. Implementation checklist for Claude", styles, "Section"))
    checklist_items = [
        "Edit TOVBrachistocrone_copia.tex, not the older reference copy.",
        "Correct W(q_0) - W(x) at source line 1266.",
        "Replace the endpoint paragraph at lines 1195-1200 with the regularized derivation.",
        "Add the explicit common-correction cancellation after Eq. (63).",
        "Break Eq. (67) with an aligned environment and recompile in two-column mode.",
        "Enable hidelinks and confirm that no colored PDF rectangles remain.",
        "Qualify every occurrence of 'any equation of state' in the abstract, Sec. X, and conclusion.",
        "Add the fixed-spatial-endpoints clarification in Sec. V.",
        "Revise the Fig. 7 validation claim and the Fig. 1 ordering language.",
        "Standardize mu in code and generated labels to r_s/R = 2M/R for paper-facing outputs.",
        "Recompile, render every page, inspect Eq. (67), references, captions, and page breaks.",
        "Run the numerical scripts from tag v1.0-submission and compare regenerated figure hashes or data tables with the submitted assets.",
    ]
    for index, checklist_item in enumerate(checklist_items, start=1):
        story.append(paragraph(f"<b>{index}.</b>&nbsp; {checklist_item}", styles, "Checklist"))

    story.append(Spacer(1, 6 * mm))
    story.append(
        callout(
            "Submission readiness after these changes",
            "Once the P0 items are fixed, the manuscript's main analytical claim is internally consistent and the "
            "constant-density benchmark is formally correct. The remaining P1 and P2 changes mainly improve the "
            "precision of the universality claim, anticipate referee questions, and make the paper easier to navigate.",
            styles,
            "ok",
        )
    )
    return story


def main() -> None:
    register_fonts()
    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
    styles = build_styles()
    document = BaseDocTemplate(
        OUTPUT_PATH,
        pagesize=A4,
        leftMargin=18 * mm,
        rightMargin=18 * mm,
        topMargin=17 * mm,
        bottomMargin=19 * mm,
        title="TOV Brachistochrone Editorial Review and Proposed Changes",
        author="OpenAI Codex",
        subject="Editorial and mathematical review of the first-order t-tau splitting derivation",
    )
    content_frame = Frame(
        document.leftMargin,
        document.bottomMargin,
        document.width,
        document.height,
        id="content",
    )
    document.addPageTemplates(
        [PageTemplate(id="review", frames=[content_frame], onPage=page_footer)]
    )
    document.build(build_story(styles))
    print(OUTPUT_PATH)


if __name__ == "__main__":
    main()
