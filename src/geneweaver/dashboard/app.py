
from textual.app import App, ComposeResult
from textual.containers import Container, Horizontal, Vertical
from textual.widgets import (
    Footer,
    Header,
    Label,
    ProgressBar,
    Static,
)


class GeneWeaverDashboard(App):
    """GeneWeaver terminal monitoring dashboard."""

    TITLE = "GeneWeaver"
    SUB_TITLE = "GPU-Accelerated CRISPR Alignment Engine"

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
                    "GENEWEAVER SYSTEM STATUS",
                    classes="panel-title",
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

    def action_refresh(self) -> None:
        self.query_one("#status", Label).update(
            "● Dashboard refreshed"
        )


if __name__ == "__main__":
    GeneWeaverDashboard().run()
