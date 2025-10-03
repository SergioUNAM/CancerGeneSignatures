from app.ui.components import Highlight, render_highlight_pills


def test_highlight_reexport_and_render_does_not_error(monkeypatch):
    # Monkeypatch Streamlit markdown to capture inputs without rendering.
    calls = []

    def fake_markdown(content, *, unsafe_allow_html=False):
        calls.append((content, unsafe_allow_html))

    monkeypatch.setattr("streamlit.markdown", fake_markdown)

    highlight = Highlight(label="Genes válidos", value=12)
    render_highlight_pills([highlight], key="test")

    assert calls, "Expected markdown to be called when rendering highlight pills"
    html = calls[0][0]
    assert "Genes válidos" in html
    assert "12" in html
