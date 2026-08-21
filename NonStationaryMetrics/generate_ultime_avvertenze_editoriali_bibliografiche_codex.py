#!/usr/bin/env python3
"""Generate the final editorial and bibliographic handoff report for Claude."""

from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate,
    Flowable,
    Frame,
    HRFlowable,
    KeepTogether,
    PageBreak,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
)


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "output" / "pdf" / "CQG_ultime_avvertenze_editoriali_bibliografiche_codex.pdf"

NAVY = colors.HexColor("#15253F")
BLUE = colors.HexColor("#2D6A9F")
TEAL = colors.HexColor("#247A78")
RED = colors.HexColor("#A53D3D")
AMBER = colors.HexColor("#B77812")
INK = colors.HexColor("#202631")
MUTED = colors.HexColor("#5F6875")
PALE_BLUE = colors.HexColor("#EAF2F8")
PALE_RED = colors.HexColor("#FBECEC")
PALE_AMBER = colors.HexColor("#FFF4DC")
PALE_GREEN = colors.HexColor("#EAF6F1")
LIGHT = colors.HexColor("#F5F7FA")
LINE = colors.HexColor("#D7DEE7")


def register_fonts():
    base = Path("/System/Library/Fonts/Supplemental")
    fonts = {
        "Arial": base / "Arial.ttf",
        "Arial-Bold": base / "Arial Bold.ttf",
        "Arial-Italic": base / "Arial Italic.ttf",
        "Arial-BoldItalic": base / "Arial Bold Italic.ttf",
        "Andale": base / "Andale Mono.ttf",
    }
    for name, path in fonts.items():
        if path.exists():
            pdfmetrics.registerFont(TTFont(name, str(path)))
    pdfmetrics.registerFontFamily(
        "Arial",
        normal="Arial",
        bold="Arial-Bold",
        italic="Arial-Italic",
        boldItalic="Arial-BoldItalic",
    )


class AccentBar(Flowable):
    def __init__(self, width, height=3 * mm, color=BLUE):
        super().__init__()
        self.width = width
        self.height = height
        self.color = color

    def draw(self):
        self.canv.setFillColor(self.color)
        self.canv.roundRect(0, 0, self.width, self.height, 1.4 * mm, fill=1, stroke=0)


class NumberedCanvasMixin:
    pass


def page_frame(canvas, doc):
    canvas.saveState()
    canvas.setTitle("Ultime avvertenze editoriali e bibliografiche - Paper I")
    canvas.setAuthor("Codex - editorial review")
    canvas.setSubject("Handoff operativo per la revisione CQG di paper1-12")
    canvas.setKeywords("CQG, Paper I, editorial review, bibliography, Claude handoff")
    w, h = A4
    if doc.page > 1:
        canvas.setStrokeColor(LINE)
        canvas.setLineWidth(0.5)
        canvas.line(19 * mm, h - 17 * mm, w - 19 * mm, h - 17 * mm)
        canvas.setFont("Arial", 8.2)
        canvas.setFillColor(MUTED)
        canvas.drawString(19 * mm, h - 13 * mm, "PAPER I - CQG FINAL HANDOFF")
        canvas.drawRightString(w - 19 * mm, h - 13 * mm, "paper1-12")
    canvas.setStrokeColor(LINE)
    canvas.setLineWidth(0.5)
    canvas.line(19 * mm, 15 * mm, w - 19 * mm, 15 * mm)
    canvas.setFont("Arial", 8.2)
    canvas.setFillColor(MUTED)
    canvas.drawString(19 * mm, 10 * mm, "Codex editorial review - 3 August 2026")
    canvas.drawRightString(w - 19 * mm, 10 * mm, f"Page {doc.page}")
    canvas.restoreState()


def build_styles():
    styles = getSampleStyleSheet()
    styles.add(
        ParagraphStyle(
            name="CoverKicker",
            fontName="Arial-Bold",
            fontSize=10,
            leading=13,
            textColor=BLUE,
            spaceAfter=4 * mm,
            tracking=0.8,
        )
    )
    styles.add(
        ParagraphStyle(
            name="CoverTitle",
            fontName="Arial-Bold",
            fontSize=25,
            leading=29,
            textColor=NAVY,
            spaceAfter=4 * mm,
        )
    )
    styles.add(
        ParagraphStyle(
            name="CoverSub",
            fontName="Arial",
            fontSize=12,
            leading=17,
            textColor=MUTED,
            spaceAfter=7 * mm,
        )
    )
    styles.add(
        ParagraphStyle(
            name="H1x",
            fontName="Arial-Bold",
            fontSize=17,
            leading=21,
            textColor=NAVY,
            spaceBefore=3 * mm,
            spaceAfter=3 * mm,
            keepWithNext=True,
        )
    )
    styles.add(
        ParagraphStyle(
            name="H2x",
            fontName="Arial-Bold",
            fontSize=12.2,
            leading=15,
            textColor=BLUE,
            spaceBefore=3.5 * mm,
            spaceAfter=1.8 * mm,
            keepWithNext=True,
        )
    )
    styles.add(
        ParagraphStyle(
            name="Bodyx",
            fontName="Arial",
            fontSize=9.6,
            leading=14.1,
            textColor=INK,
            alignment=TA_LEFT,
            spaceAfter=2.2 * mm,
        )
    )
    styles.add(
        ParagraphStyle(
            name="Smallx",
            fontName="Arial",
            fontSize=8.3,
            leading=11.5,
            textColor=MUTED,
            spaceAfter=1.5 * mm,
        )
    )
    styles.add(
        ParagraphStyle(
            name="Bulletx",
            fontName="Arial",
            fontSize=9.4,
            leading=13.6,
            leftIndent=5 * mm,
            firstLineIndent=-3.5 * mm,
            bulletIndent=1 * mm,
            textColor=INK,
            spaceAfter=1.8 * mm,
        )
    )
    styles.add(
        ParagraphStyle(
            name="CalloutTitle",
            fontName="Arial-Bold",
            fontSize=10.2,
            leading=13,
            textColor=NAVY,
            spaceAfter=1.3 * mm,
        )
    )
    styles.add(
        ParagraphStyle(
            name="CalloutBody",
            fontName="Arial",
            fontSize=9.2,
            leading=13.2,
            textColor=INK,
        )
    )
    styles.add(
        ParagraphStyle(
            name="Quote",
            fontName="Arial-Italic",
            fontSize=9.1,
            leading=13.2,
            leftIndent=5 * mm,
            rightIndent=4 * mm,
            textColor=INK,
            borderColor=LINE,
            borderWidth=0.6,
            borderPadding=3 * mm,
            backColor=LIGHT,
            spaceBefore=1.5 * mm,
            spaceAfter=2.5 * mm,
        )
    )
    styles.add(
        ParagraphStyle(
            name="CodeBlockx",
            fontName="Andale",
            fontSize=7.8,
            leading=11.2,
            leftIndent=3.5 * mm,
            rightIndent=3.5 * mm,
            borderColor=LINE,
            borderWidth=0.6,
            borderPadding=3 * mm,
            backColor=colors.white,
            textColor=INK,
            spaceBefore=1.5 * mm,
            spaceAfter=2.5 * mm,
        )
    )
    styles.add(
        ParagraphStyle(
            name="TableHead",
            fontName="Arial-Bold",
            fontSize=8.5,
            leading=10.5,
            textColor=colors.white,
            alignment=TA_LEFT,
        )
    )
    styles.add(
        ParagraphStyle(
            name="TableCell",
            fontName="Arial",
            fontSize=8.3,
            leading=11.3,
            textColor=INK,
        )
    )
    return styles


def P(text, style):
    return Paragraph(text, style)


def bullet(text, styles):
    return Paragraph(f"<bullet>&#8226;</bullet>{text}", styles["Bulletx"])


def callout(title, body, styles, bg=PALE_BLUE, accent=BLUE):
    inner = Table(
        [[P(title, styles["CalloutTitle"])], [P(body, styles["CalloutBody"])]],
        colWidths=[160 * mm],
        hAlign="LEFT",
    )
    inner.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), bg),
                ("BOX", (0, 0), (-1, -1), 0.7, accent),
                ("LINEBEFORE", (0, 0), (0, -1), 4, accent),
                ("LEFTPADDING", (0, 0), (-1, -1), 8),
                ("RIGHTPADDING", (0, 0), (-1, -1), 8),
                ("TOPPADDING", (0, 0), (-1, 0), 7),
                ("BOTTOMPADDING", (0, -1), (-1, -1), 8),
            ]
        )
    )
    return inner


def issue(title, location, problem, action, styles, priority="MUST FIX"):
    color = RED if priority == "MUST FIX" else AMBER
    bg = PALE_RED if priority == "MUST FIX" else PALE_AMBER
    rows = [
        [P(priority, styles["TableHead"]), P(title, styles["CalloutTitle"])],
        [P("LOCATION", styles["TableHead"]), P(location, styles["TableCell"])],
        [P("PROBLEM", styles["TableHead"]), P(problem, styles["TableCell"])],
        [P("ACTION", styles["TableHead"]), P(action, styles["TableCell"])],
    ]
    table = Table(rows, colWidths=[26 * mm, 134 * mm], hAlign="LEFT", splitByRow=1)
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (0, -1), color),
                ("BACKGROUND", (1, 0), (1, 0), bg),
                ("BACKGROUND", (1, 1), (1, -1), colors.white),
                ("TEXTCOLOR", (0, 0), (0, -1), colors.white),
                ("BOX", (0, 0), (-1, -1), 0.7, color),
                ("INNERGRID", (0, 1), (-1, -1), 0.35, LINE),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("TOPPADDING", (0, 0), (-1, -1), 6),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
            ]
        )
    )
    return KeepTogether([table, Spacer(1, 3.5 * mm)])


def checklist_row(status, item, styles):
    return [P(status, styles["TableCell"]), P(item, styles["TableCell"])]


def make_story(styles):
    story = []
    usable = 160 * mm
    story.extend(
        [
            Spacer(1, 12 * mm),
            AccentBar(48 * mm, color=BLUE),
            Spacer(1, 7 * mm),
            P("FINAL CQG HANDOFF", styles["CoverKicker"]),
            P("Ultime avvertenze editoriali e bibliografiche", styles["CoverTitle"]),
            P(
                "Istruzioni operative per Claude sulla build <b>paper1-12</b>. "
                "Il nucleo scientifico e la correzione adiabatica completa sono ora sostanzialmente "
                "presentabili; questo documento isola gli ultimi interventi necessari prima della submission.",
                styles["CoverSub"],
            ),
            callout(
                "VERDETTO",
                "Non inviare ancora la build attuale. La revisione non richiede piu una ristrutturazione "
                "scientifica generale, ma contiene tre blocchi reali: una contraddizione sulle degenerazioni, "
                "un conflitto semantico fra J_c e J_deg, e una tabella sovrapposta. Dopo questi interventi e "
                "una breve pulizia Paper I/Paper II, la versione puo essere considerata pronta per il referee.",
                styles,
                bg=PALE_RED,
                accent=RED,
            ),
            Spacer(1, 6 * mm),
        ]
    )

    summary_data = [
        [P("AREA", styles["TableHead"]), P("STATUS", styles["TableHead"]), P("DECISION", styles["TableHead"])],
        [P("Core scientifico", styles["TableCell"]), P("Substantially ready", styles["TableCell"]), P("Preservare Teorema I.5 e Appendix C", styles["TableCell"])],
        [P("Notazione", styles["TableCell"]), P("Blocking", styles["TableCell"]), P("Separare J_cap, J_deg e J_sep", styles["TableCell"])],
        [P("Bibliografia", styles["TableCell"]), P("Good, final cleanup", styles["TableCell"]), P("43/43 citate; completare metadati e DOI", styles["TableCell"])],
        [P("Produzione PDF", styles["TableCell"]), P("Blocking", styles["TableCell"]), P("Riparare Table A2 e i float", styles["TableCell"])],
    ]
    summary = Table(summary_data, colWidths=[39 * mm, 38 * mm, 83 * mm], repeatRows=1, hAlign="LEFT")
    summary.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), NAVY),
                ("GRID", (0, 0), (-1, -1), 0.5, LINE),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, LIGHT]),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("TOPPADDING", (0, 0), (-1, -1), 6),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
            ]
        )
    )
    story.extend([summary, Spacer(1, 5 * mm)])
    story.append(
        P(
            "<b>Documento ispezionato:</b> 1-paper1-12.pdf, 37 pagine, generato il 3 August 2026. "
            "Controllo eseguito sul testo estratto e sul rendering integrale pagina per pagina.",
            styles["Smallx"],
        )
    )
    story.append(PageBreak())

    story.extend([P("1. Correzioni obbligatorie", styles["H1x"]), AccentBar(usable, color=RED), Spacer(1, 4 * mm)])
    story.append(
        issue(
            "Contraddizione fra freezing, turning point e separatrix",
            "PDF pages 4 and 9; section 2.1 and the domain discussion after Proposition I.3.",
            "Page 4 identifies the freezing surface with the turning-point/separatrix locus. Page 9 "
            "states, correctly, that freezing, turning and the double-root separatrix are distinct degeneracies. "
            "Both statements cannot remain.",
            "Keep the distinction stated on page 9. Replace the parenthesis on page 4 with: "
            "<i>'At the freezing surface the entire indicatrix collapses, R -> 0, and strict convexity is lost. "
            "This degeneration is distinct from a radial turning point and from a spectral separatrix, which are "
            "treated separately below.'</i>",
            styles,
        )
    )
    story.append(
        issue(
            "J_c and J_deg currently denote different objects",
            "Section 4.2 (page 14), section 4.4 (pages 17-20), Appendix B (pages 27-31), and Appendix C (pages 34-35).",
            "In the main text J_c(v_0) is the dynamical penetration/capture threshold. J_deg is the algebraic "
            "value at which the sextic develops a negative-radius double root. Appendix B then writes "
            "J = J_c = J_deg, and Appendix C calls J -> J_c a separatrix limit. This reuses one symbol for "
            "physically inequivalent loci.",
            "Do not perform a blind global replacement. Use J_cap(v_0) or J_pen(v_0) for the physical Vaidya "
            "capture threshold; use J_deg for the algebraic genus degeneration at r_d < 0; reserve J_sep for a "
            "physical accessible double-root separatrix (primarily Paper II). In Vaidya appendices, write "
            "'algebraic degeneration locus' unless a physical separatrix has actually been established.",
            styles,
        )
    )
    story.append(
        issue(
            "Table A2 contains overlapping text",
            "PDF page 26, last row of Table A2.",
            "The script name vaidya_first_order_offshell.py overlaps the residual entry 'slope 2.00'. The row is "
            "visibly broken at normal reading scale.",
            "Rebuild the table with a wrapping script column, a smaller monospaced size, or a two-line final row. "
            "Render the page and confirm that the filename, residual and SHA-256 line do not collide.",
            styles,
        )
    )
    story.append(
        issue(
            "Paper II material still interrupts the Vaidya argument",
            "Table 1 and Figure 1; pages 4, 9, 13-14, 17-21; repeated appendix forward references.",
            "Paper I still carries a Thakurta-Kerr notation block, a full Thakurta-Kerr indicatrix panel, an "
            "ergosphere domain excursus and many implementation-level comparisons. The sentence on page 21 "
            "equating the Vaidya horizon-weight sector with an ergosphere-penetration window is especially "
            "misleading because spherical Vaidya has no ergosphere.",
            "Retain only short forward references to Paper II. Delete the Thakurta-Kerr block from Table 1, remove "
            "or relocate the third panel of Figure 1, compress the ergosphere discussion to one domain caveat, and "
            "replace the page-21 analogy with: <i>'Thus the advanced/retarded frozen-clock asymmetry isolates the "
            "contribution of the horizon-anchored weight-two sector.'</i>",
            styles,
            priority="SHOULD FIX",
        )
    )

    story.extend([P("2. Bibliografia: stato e interventi", styles["H1x"]), AccentBar(usable, color=TEAL), Spacer(1, 4 * mm)])
    story.append(
        callout(
            "STATO POSITIVO DA PRESERVARE",
            "The current bibliography contains 43 entries and all 43 are explicitly cited in the manuscript. "
            "The former nocite-style inflation is no longer present. Do not reintroduce uncited entries, and do not "
            "remove references [8]-[11], which now supply the missing relativistic brachistochrone/Fermat context.",
            styles,
            bg=PALE_GREEN,
            accent=TEAL,
        )
    )
    story.append(Spacer(1, 4 * mm))

    bib_rows = [
        [P("ITEM", styles["TableHead"]), P("WHAT CLAUDE SHOULD DO", styles["TableHead"])],
        [P("References heading", styles["TableCell"]), P("Insert a visible <b>References</b> heading before item [1]. At present [1] follows the ORCID line without a section heading.", styles["TableCell"])],
        [P("DOI pass", styles["TableCell"]), P("Where a DOI exists, add it after verifying it against the publisher or Crossref. IOP asks for a permanent identifier when available. Never infer a DOI from a title.", styles["TableCell"])],
        [P("Reference [11]", styles["TableCell"]), P("Keep as a preprint unless a later published record is found. Verified title: <i>On the formulations of the Fermat principle in general relativity and beyond</i>; authors Erasmo Caponio and Miguel Angel Javaloyes; arXiv:2605.01532.", styles["TableCell"])],
        [P("Reference [15]", styles["TableCell"]), P("Expand the incomplete Memoirs citation. Official metadata: Caponio E, Javaloyes M A and Sánchez M, <i>Wind Finslerian Structures: From Zermelo's Navigation to the Causality of Spacetimes</i>, <i>Mem. Amer. Math. Soc.</i> <b>300</b>, no. 1501 (2024), 121 pp.", styles["TableCell"])],
        [P("Reference [31]", styles["TableCell"]), P("Replace the generic 'arXiv preprint' wording with a locatable entry: Brown F C S and Levin A 2011 <i>Multiple Elliptic Polylogarithms</i> arXiv:1110.6917.", styles["TableCell"])],
        [P("Reference [12]", styles["TableCell"]), P("Protect proper-noun capitalization in BibTeX: {Thakurta}, {Kerr}, and {Paper II}. Keep the Zenodo DOI and avoid phrasing that reads like an unpublished result being used as external validation of Paper I.", styles["TableCell"])],
        [P("Reference [32]", styles["TableCell"]), P("Keep the Zenodo DOI as the primary persistent identifier. The GitHub URL and exact commit may remain as secondary reproducibility information; prevent awkward line breaks around the DOI and commit hash.", styles["TableCell"])],
        [P("Reference [39]", styles["TableCell"]), P("If the patched abelfunctions tree has an archived commit or snapshot identifier, cite that exact identifier. Typeset is_* as code so the asterisk is not interpreted as emphasis or a wildcard outside its API context.", styles["TableCell"])],
        [P("Consistency", styles["TableCell"]), P("Use one consistent convention for preprints: either 'arXiv:xxxx.xxxxx' or '(Preprint xxxx.xxxxx)'. Preserve accents in Sanchez and other author names through BibTeX braces/UTF-8.", styles["TableCell"])],
    ]
    bib_table = Table(bib_rows, colWidths=[35 * mm, 125 * mm], repeatRows=1, hAlign="LEFT")
    bib_table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), TEAL),
                ("GRID", (0, 0), (-1, -1), 0.45, LINE),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, LIGHT]),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("TOPPADDING", (0, 0), (-1, -1), 6),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
            ]
        )
    )
    story.extend([bib_table, Spacer(1, 4 * mm)])
    story.append(
        P(
            "<b>Important:</b> article titles are optional under the general IOP journal style, but journal, volume, "
            "page/article number and DOI (when available) should make each item independently locatable. Preprints "
            "must at least carry authors, year and preprint number.",
            styles["Bodyx"],
        )
    )

    story.extend([P("3. Produzione e leggibilita", styles["H1x"]), AccentBar(usable, color=AMBER), Spacer(1, 4 * mm)])
    for item in [
        "Figure 1 on page 7 is too small at normal reading scale. If the Thakurta-Kerr panel is removed, enlarge the two Paper-I panels and harmonize legend/font sizes.",
        "Pages 21, 31 and 35 are largely empty because of float/page-break behavior. Reflow the text and floats; this should remove roughly two to three pages without cutting scientific content.",
        "Figure 9 should state the epsilon fit interval, number of sampled points, residual norm and fitted slopes (with uncertainties if readily available). Preserve the current O(epsilon) versus O(epsilon^2) comparison.",
        "Fill PDF metadata: title, author, subject and keywords are currently blank. All fonts are already embedded, so preserve that successful part of the build.",
        "Run a final reference/cross-reference check after every notation change: no ?? markers, no stale equation numbers, and no accidental renumbering of Theorem I.5, equation (22), Figure 9 or Appendices B-C.",
    ]:
        story.append(bullet(item, styles))

    story.extend([P("4. Materiale da non alterare", styles["H1x"]), AccentBar(usable, color=TEAL), Spacer(1, 4 * mm)])
    for item in [
        "Keep Theorem I.5 in the main text. It is the central complete on-shell plus off-shell result and should not be buried again in Appendix C.",
        "Keep the explicit caveat that global minimization is conditional and not automatic. The current abstract and section 2.2 correctly distinguish PMP extremals, local checks and an HJB certificate.",
        "Keep the statement that the first-order correction is non-uniform at the turning point and is verified on compact sub-arcs.",
        "Keep the hierarchy of special-function classes and the explicit admission that weight-one reducibility and a canonical single-valued completion remain conjectural/open.",
        "Keep Figure 9 and the slope-2 validation. Only improve its caption and physical size.",
        "Keep the curated references [8]-[11], the Zenodo reproducibility citation [32], and the transparent AI-use disclosure unless the journal requests different placement.",
    ]:
        story.append(bullet(item, styles))

    story.extend([P("5. Checklist di accettazione", styles["H1x"]), AccentBar(usable, color=BLUE), Spacer(1, 4 * mm)])
    checklist = [
        checklist_row("[ ]", "Page 4 no longer identifies freezing with turning/separatrix.", styles),
        checklist_row("[ ]", "J_cap/J_pen, J_deg and J_sep have distinct meanings everywhere.", styles),
        checklist_row("[ ]", "Appendices B-C no longer call the negative-radius algebraic degeneration a physical separatrix.", styles),
        checklist_row("[ ]", "Table A2 has no overlapping text at 100% zoom.", styles),
        checklist_row("[ ]", "The Vaidya discussion no longer depends on an ergosphere analogy.", styles),
        checklist_row("[ ]", "A visible References heading appears before item [1].", styles),
        checklist_row("[ ]", "References [11], [15], [31], [32] and [39] have complete, locatable metadata.", styles),
        checklist_row("[ ]", "All available DOIs added only after verification; all 43 references still cited.", styles),
        checklist_row("[ ]", "Figure 1 is legible; pages 21, 31 and 35 have been reflowed.", styles),
        checklist_row("[ ]", "Figure 9 caption specifies fit domain and residual definition.", styles),
        checklist_row("[ ]", "PDF metadata populated; fonts remain embedded.", styles),
        checklist_row("[ ]", "Fresh full-PDF rendering shows no collisions, clipping, empty float pages or stale links.", styles),
    ]
    check_table = Table(checklist, colWidths=[13 * mm, 147 * mm], hAlign="LEFT")
    check_table.setStyle(
        TableStyle(
            [
                ("GRID", (0, 0), (-1, -1), 0.45, LINE),
                ("ROWBACKGROUNDS", (0, 0), (-1, -1), [colors.white, LIGHT]),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("ALIGN", (0, 0), (0, -1), "CENTER"),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("TOPPADDING", (0, 0), (-1, -1), 7),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 7),
            ]
        )
    )
    story.extend([check_table, Spacer(1, 6 * mm)])

    story.extend([P("6. Prompt pronto da incollare in Claude", styles["H1x"]), AccentBar(usable, color=NAVY), Spacer(1, 4 * mm)])
    prompt = (
        "Revise Paper I using the attached Codex handoff as an authoritative scope list. Fix every MUST FIX item "
        "and implement the bibliographic cleanup, but do not rewrite or rederive the scientific core. Preserve "
        "Theorem I.5 in the main text, equation (22), Appendix C, Figure 9, all numerical residuals, and the "
        "conditional wording on global optimality. Distinguish semantically between the physical Vaidya capture "
        "threshold J_cap(v_0), the algebraic degeneration J_deg at r_d < 0, and any genuinely physical separatrix "
        "J_sep. Do not perform blind global replacements. Remove the false identification of freezing with "
        "turning/separatrix, repair Table A2, reduce Paper-II leakage, add a References heading, complete the "
        "specified bibliography entries, and add DOIs only when independently verified. After editing, compile "
        "the PDF and inspect every page at normal zoom. Report each changed source location, any equation/figure "
        "renumbering, and the final citation audit count. Do not change any coefficient, sign, residual or special-"
        "function claim unless you identify and explain a concrete inconsistency first."
    )
    story.append(P(prompt, styles["CodeBlockx"]))

    story.extend([P("7. Fonti usate per il controllo bibliografico", styles["H1x"]), AccentBar(usable, color=TEAL), Spacer(1, 4 * mm)])
    sources = [
        ("IOP Publishing, Style guide for journal articles", "https://publishingsupport.iopscience.iop.org/questions/style-guide-journal-articles/"),
        ("Caponio, Javaloyes and Sánchez, official AMS Memoirs record", "https://bookstore.ams.org/memo-300-1501"),
        ("Brown and Levin, Multiple Elliptic Polylogarithms", "https://arxiv.org/abs/1110.6917"),
        ("Caponio and Javaloyes, Fermat principle survey", "https://arxiv.org/abs/2605.01532"),
    ]
    source_cells = []
    for label, url in sources:
        source_cells.append(P(f'<link href="{url}" color="#2D6A9F"><b>{label}</b></link><br/><font size="7.1" color="#5F6875">{url}</font>', styles["Smallx"]))
    source_table = Table(
        [[source_cells[0], source_cells[1]], [source_cells[2], source_cells[3]]],
        colWidths=[78 * mm, 78 * mm],
        hAlign="LEFT",
    )
    source_table.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (-1, -1), 8),
                ("TOPPADDING", (0, 0), (-1, -1), 2),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
            ]
        )
    )
    story.append(source_table)
    return story


def main():
    register_fonts()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    styles = build_styles()
    doc = BaseDocTemplate(
        str(OUTPUT),
        pagesize=A4,
        leftMargin=19 * mm,
        rightMargin=19 * mm,
        topMargin=22 * mm,
        bottomMargin=20 * mm,
        title="Ultime avvertenze editoriali e bibliografiche - Paper I",
        author="Codex",
        subject="Handoff operativo per Claude sulla build paper1-12",
    )
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id="main")
    doc.addPageTemplates([PageTemplate(id="all", frames=frame, onPage=page_frame)])
    doc.build(make_story(styles))
    print(OUTPUT)


if __name__ == "__main__":
    main()
