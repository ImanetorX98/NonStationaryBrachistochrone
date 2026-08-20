#!/usr/bin/env bash
# Rebuild Paper I from a clean checkout.
#
# Figures are not tracked (the repository ignores *.pdf and *.png), so they must
# be regenerated from their scripts before LaTeX can find them. This script does
# that and then compiles, so that a fresh clone is turnkey.
#
#   ./build_paper1.sh            # figures + manuscript
#   ./build_paper1.sh --figures  # figures only
#   ./build_paper1.sh --paper    # manuscript only (figures assumed present)
#
# Requires: python3 with numpy/scipy/sympy/mpmath/matplotlib, and a TeX
# distribution with latexmk. The symbolic verification suite additionally wants
# wolframscript; see verification/README.md.

set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
MODE="${1:-all}"

run_figs () {
  echo "== regenerating figures into paper/Immagini =="
  local failed=0
  # one entry per generating script; several produce more than one figure
  local scripts=(
    "ThakurtaMetric/genera_figure_thakurta.py"        # fig_indicatrici
    "FLRWmetric/genera_figure_flrw.py"                # fig_flrw_cinematica, _variazionale
    "FLRWmetric/genera_fig_flrw_curvatura.py"         # fig_flrw_curvatura
    "VaidyaMetric/kodama_conservation.py"             # fig_kodama_conservazione
    "VaidyaMetric/genera_figure_vaidya.py"            # fig_vaidya_traiettorie, _kerr_a0
    "VaidyaMetric/plunge_vaidya_t_tau.py"             # fig_vaidya_plunge_t_tau
    "VaidyaMetric/inversione_fisica.py"               # fig_vaidya_no_inversione_evaporazione
    "VaidyaMetric/vaidya_brachistochrone_vparam.py"   # fig_vaidya_bounce, fig_vaidya_timing
    "VaidyaMetric/vaidya_first_order_offshell_plot_codex.py"  # fig_vaidya_offshell
    "VaidyaMetric/verifica_minimo_brachi.py"          # fig_verifica_minimo_brachi
  )
  for s in "${scripts[@]}"; do
    printf '  %-56s ' "$s"
    if [ ! -f "$ROOT/$s" ]; then echo "MISSING"; failed=$((failed+1)); continue; fi
    if ( cd "$ROOT/$(dirname "$s")" && python3 "$(basename "$s")" ) >/dev/null 2>&1; then
      echo "ok"
    else
      echo "FAILED"; failed=$((failed+1))
    fi
  done
  echo "  ($failed script(s) failed)"
}

collect_figs () {
  # each generator writes into its own directory; the manuscript looks in
  # paper/Immagini, so gather them there
  echo "== collecting figures into paper/Immagini =="
  mkdir -p "$ROOT/paper/Immagini"
  local n=0
  for d in "$ROOT/PaperFigures" "$ROOT/FLRWmetric/FLRWfigures" \
           "$ROOT/VaidyaMetric/Vaidyafigures" "$ROOT/ThakurtaMetric/Thakurtafigures"; do
    [ -d "$d" ] || continue
    for f in "$d"/fig_*.pdf; do
      [ -e "$f" ] || continue
      cp -f "$f" "$ROOT/paper/Immagini/" && n=$((n+1))
    done
  done
  echo "  copied $n file(s)"
}

check_figs () {
  local missing=0
  for f in fig_indicatrici fig_flrw_cinematica fig_flrw_variazionale \
           fig_flrw_curvatura fig_kodama_conservazione fig_vaidya_traiettorie \
           fig_vaidya_plunge_t_tau fig_vaidya_no_inversione_evaporazione \
           fig_vaidya_bounce fig_vaidya_timing fig_vaidya_offshell \
           fig_vaidya_kerr_a0 fig_verifica_minimo_brachi; do
    [ -f "$ROOT/paper/Immagini/$f.pdf" ] || { echo "  missing: $f.pdf"; missing=$((missing+1)); }
  done
  echo "== $((13-missing))/13 figures present =="
  return 0
}

run_paper () {
  echo "== compiling paper1_JMP.tex =="
  ( cd "$HERE" && latexmk -pdf -interaction=nonstopmode -halt-on-error paper1_JMP.tex >/dev/null 2>&1 )
  if [ -f "$HERE/paper1_JMP.pdf" ]; then
    echo "  built paper1_JMP.pdf"
  else
    echo "  BUILD FAILED - see paper1_JMP.log"; return 1
  fi
}

case "$MODE" in
  --figures) run_figs; collect_figs; check_figs ;;
  --paper)   collect_figs; check_figs; run_paper ;;
  *)         run_figs; collect_figs; check_figs; run_paper ;;
esac
