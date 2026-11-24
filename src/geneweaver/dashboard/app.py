
from textual.app import App, ComposeResult
from textual.containers import Container, Horizontal, Vertical
from textual.widgets import (
    Footer,
    Header,
    Label,
    ProgressBar,
    Static,
    Select,
)
from pathlib import Path
from geneweaver.fasta import iter_fasta_chunks


class GeneWeaverDashboard(App):
    """GeneWeaver terminal monitoring dashboard."""

    TITLE = "GeneWeaver"
    SUB_TITLE = "GPU-Accelerated CRISPR Alignment Engine"
    DATASETS = [
        ("Mock Genome", "data/mock_genome.fasta"),
        ("100 KB Genome", "data/benchmarks/genome_100kb.fasta"),
        ("1 MB Genome", "data/benchmarks/genome_1mb.fasta"),
        ("10 MB Genome", "data/benchmarks/genome_10mb.fasta"),
    ]

    CSS = """
    Screen {
        background: #101827;
    }

    #main-container {
        padding: 1 2;
        height: 1fr;
    }

    .panel {
        border: round #2563eb;
        background: #172337;
        padding: 1 2;
        margin-bottom: 1;
        height: auto;
    }

    .panel-title {
        color: #60a5fa;
        text-style: bold;
        margin-bottom: 1;
    }

    .metric {
        color: #e2e8f0;
        margin-bottom: 1;
    }

    #status {
        color: #4ade80;
        text-style: bold;
    }

    #progress {
        margin-top: 1;
        margin-bottom: 1;
    }

    #footer-note {
        color: #94a3b8;
        margin-top: 1;
    }
    """

    BINDINGS = [
        ("q", "quit", "Quit"),
        ("r", "refresh", "Refresh"),
    ]

    def compose(self) -> ComposeResult:
        yield Header()

        with Container(id="main-container"):

            with Vertical(classes="panel"):
                yield Label(
    "Genome size: Waiting for input",
    id="genome-size",
    classes="metric",
)
                yield Select(
        self.DATASETS,
        prompt="Choose a genome",
        id="dataset-select",
    )
                yield Label(
                    "● System initialized",
                    id="status",
                )
                yield Label(
                    "Processing engine: Ready",
                    classes="metric",
                )
                yield Label(
                    "GPU backend: OpenCL",
                    classes="metric",
                )

            with Horizontal():

                with Vertical(classes="panel"):
                    yield Label(
                        "GENOME INFORMATION",
                        classes="panel-title",
                    )
                    yield Label(
                        "Dataset: Not selected",
                        id="dataset",
                        classes="metric",
                    )
                    yield Label(
                        "Genome size: Waiting for input",
                        classes="metric",
                    )
                    yield Label(
                        "Sequence type: DNA",
                        classes="metric",
                    )

                with Vertical(classes="panel"):
                    yield Label(
                        "ALIGNMENT METRICS",
                        classes="panel-title",
                    )
                    yield Label(
                        "CPU alignment: Ready",
                        classes="metric",
                    )
                    yield Label(
                        "GPU alignment: Ready",
                        classes="metric",
                    )
                    yield Label(
                        "Off-target hits: 0",
                        classes="metric",
                    )

            with Vertical(classes="panel"):
                yield Label(
                    "PROCESSING PROGRESS",
                    classes="panel-title",
                )
                yield ProgressBar(
                    total=100,
                    show_eta=False,
                    id="progress",
                )
                yield Label(
                    "Waiting for alignment workload...",
                    id="progress-label",
                    classes="metric",
                )

            yield Static(
                "GeneWeaver | Week 2 Day 6 | Dashboard Foundation",
                id="footer-note",
            )

            yield Footer()

    def on_select_changed(self, event: Select.Changed) -> None:
        """Load selected FASTA metadata into the dashboard."""

        if event.value is Select.BLANK:
            return

        selected_path = Path(str(event.value))

        if not selected_path.is_file():
            self.query_one("#status", Label).update(
                "File error: FASTA dataset not found"
            )
            self.query_one("#genome-size", Label).update(
                "Genome size: Unavailable"
            )
            return

        selected_name = next(
            (
                name
                for name, path in self.DATASETS
                if path == str(selected_path)
            ),
            selected_path.name,
        )

        self.query_one("#dataset", Label).update(
            f"Dataset: {selected_name}"
        )

        genome_size = sum(
            chunk.core_length
            for chunk in iter_fasta_chunks(selected_path)
        )

        self.query_one("#genome-size", Label).update(
            f"Genome size: {genome_size:,} bases"
        )

        self.query_one("#status", Label).update(
            "● Dataset loaded successfully"
        )


if __name__ == "__main__":
    GeneWeaverDashboard().run()
