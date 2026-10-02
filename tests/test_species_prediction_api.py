from pathlib import Path


def test_cli_uses_current_species_prediction_api():
    cli = (Path(__file__).parents[1] / "scripts" / "malaria-profiler").read_text()

    assert "species_prediction.species" not in cli
    assert "species_prediction.taxa" in cli
