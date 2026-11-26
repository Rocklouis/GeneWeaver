
from textual.app import App, ComposeResult
from textual import work
from textual.containers import Horizontal, Vertical, VerticalScroll
from textual.widgets import (
    Footer,
    Header,
    Label,
    ProgressBar,
    Static,
    Select,
    Input,
    Button,
)
from pathlib import Path
from geneweaver.fasta import iter_fasta_chunks
from geneweaver.dashboard.processor import process_fasta_alignment


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

        with VerticalScroll(id="main-container"):

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

            with Vertical(classes="panel"):
                yield Label(
                    "ALIGNMENT CONFIGURATION",
                    classes="panel-title",
                )
                yield Label("Guide RNA sequence:", classes="metric")
                yield Input(
                    placeholder="Enter guide sequence (e.g. ACGTACGT)",
                    id="guide-input",
                )
                yield Label("Processing engine:", classes="metric")
                yield Select(
                    [
                        ("CPU", "cpu"),
                        ("GPU (OpenCL)", "gpu"),
                    ],
                    value="cpu",
                    id="engine-select",
                )
                yield Button(
                    "Start Alignment",
                    id="start-alignment",
                    variant="success",
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
                        id="genome-size-detail",
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
                        id="hits",
                        classes="metric",
                    )
                    yield Label(
                        "Chunks processed: 0",
                        id="chunks-processed",
                        classes="metric",
                    )
                    yield Label(
                        "Bases processed: 0",
                        id="bases-processed",
                        classes="metric",
                    )
                    yield Label(
                        "Alignment time: 0.00 seconds",
                        id="alignment-time",
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

        # Ignore changes from the CPU/GPU dropdown.
        if event.control.id != "dataset-select":
            return

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
            self.query_one("#genome-size-detail", Label).update(
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

        self.query_one("#genome-size-detail", Label).update(
            f"Genome size: {genome_size:,} bases"
        )

        self.query_one("#status", Label).update(
            "Dataset loaded successfully"
        )



    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Handle the Start Alignment button."""

        if event.button.id != "start-alignment":
            return

        guide = self.query_one("#guide-input", Input).value.strip()
        engine = self.query_one("#engine-select", Select).value
        dataset = self.query_one("#dataset-select", Select).value

        if not guide:
            self.query_one("#status", Label).update(
                "Please enter a guide RNA sequence."
            )
            return

        if dataset is Select.BLANK:
            self.query_one("#status", Label).update(
                "Please select a FASTA dataset."
            )
            return

        dataset_path = Path(str(dataset))

        if not dataset_path.is_file():
            self.query_one("#status", Label).update(
                "Selected FASTA dataset was not found."
            )
            return

        if engine is Select.BLANK:
            engine = "cpu"

        self.query_one("#status", Label).update(
            "Alignment running..."
        )
        self.query_one("#start-alignment", Button).disabled = True

        self.query_one("#progress", ProgressBar).update(
            total=100,
            progress=0,
        )

        self.query_one("#progress-label", Label).update(
            "Preparing alignment..."
        )

        self.run_alignment(
            str(dataset_path),
            guide,
            str(engine),
        )

    @work(thread=True)
    def run_alignment(
        self,
        fasta_path: str,
        guide: str,
        engine: str,
    ) -> None:
        """Run alignment in a background thread."""

        try:
            result = process_fasta_alignment(
                fasta_path,
                guide,
                engine=engine,
                progress_callback=lambda *values: self.call_from_thread(
                    self.update_progress,
                    *values,
                ),
            )

            self.call_from_thread(
                self.finish_alignment,
                result,
            )

        except Exception as exc:
            self.call_from_thread(
                self.fail_alignment,
                str(exc),
            )

    def update_progress(
        self,
        processed_chunks: int,
        total_chunks: int,
        processed_bases: int,
        total_bases: int,
        hits: int,
    ) -> None:
        """Update dashboard progress safely on the UI thread."""

        percentage = (
            processed_chunks / total_chunks * 100
            if total_chunks
            else 0
        )

        self.query_one("#progress", ProgressBar).update(
            total=max(total_chunks, 1),
            progress=processed_chunks,
        )

        self.query_one("#progress-label", Label).update(
            f"Processing: {percentage:.1f}%"
        )

        self.query_one("#chunks-processed", Label).update(
            f"Chunks processed: {processed_chunks}/{total_chunks}"
        )

        self.query_one("#bases-processed", Label).update(
            f"Bases processed: {processed_bases:,}/{total_bases:,}"
        )

        self.query_one("#hits", Label).update(
            f"Off-target hits: {hits}"
        )

    def finish_alignment(self, result: dict) -> None:
        """Display final alignment results."""

        chunks = result["chunks_processed"]

        self.query_one("#progress", ProgressBar).update(
            total=max(chunks, 1),
            progress=chunks,
        )

        self.query_one("#progress-label", Label).update(
            "Alignment completed successfully."
        )

        self.query_one("#chunks-processed", Label).update(
            f"Chunks processed: {chunks}"
        )

        self.query_one("#bases-processed", Label).update(
            f"Bases processed: {result['bases_processed']:,}"
        )

        self.query_one("#hits", Label).update(
            f"Off-target hits: {result['total_hits']}"
        )

        self.query_one("#alignment-time", Label).update(
            f"Alignment time: {result['elapsed_seconds']:.4f} seconds"
        )

        self.query_one("#status", Label).update(
            f"Alignment complete ({result['engine'].upper()})"
        )

        self.query_one("#start-alignment", Button).disabled = False

    def fail_alignment(self, error: str) -> None:
        """Display an alignment error."""

        self.query_one("#status", Label).update(
            f"Alignment failed: {error}"
        )

        self.query_one("#progress-label", Label).update(
            "Alignment failed."
        )

        self.query_one("#start-alignment", Button).disabled = False



if __name__ == "__main__":
    GeneWeaverDashboard().run()
