"""
app.py
Flask web application for AI4VPrimer Amplicon-Sanger Suite.
Provides a modern no-code browser interface for running Sanger primer designs.
"""

from flask import Flask, render_template, request, jsonify
import webbrowser
import threading
import os
from sanger_designer.pipeline import SangerAmpliconPipeline

app = Flask(__name__, template_folder="templates")


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/api/run", methods=["POST"])
def run_pipeline_api():
    try:
        data = request.get_json() or {}
        fasta_path = data.get("fasta_path")
        if not fasta_path or not os.path.exists(fasta_path):
            return jsonify({
                "status": "ERROR",
                "message": f"Input FASTA file not found: {fasta_path}"
            }), 400

        fwd_pcr = data.get("fwd_pcr")
        rev_pcr = data.get("rev_pcr")
        target_start = data.get("target_start")
        target_end = data.get("target_end")
        min_cov = float(data.get("min_cov", 80.0))
        min_tm = float(data.get("min_tm", 55.0))
        max_tm = float(data.get("max_tm", 60.0))
        out_path = data.get("out_path")

        pipeline = SangerAmpliconPipeline(
            fasta_path=fasta_path,
            fwd_pcr_primer=fwd_pcr,
            rev_pcr_primer=rev_pcr,
            target_subregion_start=target_start,
            target_subregion_end=target_end,
            min_coverage_pct=min_cov,
            min_tm=min_tm,
            max_tm=max_tm,
            output_report_path=out_path
        )

        res = pipeline.run()

        return jsonify({
            "status": "SUCCESS",
            "report": res["report_content"],
            "pcr_info": res["pcr_info"],
            "tiling_info": res["tiling_info"]
        })

    except Exception as e:
        return jsonify({
            "status": "ERROR",
            "message": str(e)
        }), 500


def open_browser():
    webbrowser.open_new("http://127.0.0.1:5001")


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5001))
    print(f"Starting Amplicon-Sanger Suite on http://127.0.0.1:{port}")
    # Run Flask
    app.run(host="127.0.0.1", port=port, debug=False)
