import simplesam


def test_version_import():
    """Test that simplesam can be imported and has a version."""
    # This test ensures that the version retrieval works without pkg_resources errors
    assert hasattr(simplesam, '__version__')
    assert isinstance(simplesam.__version__, str)
    assert len(simplesam.__version__) > 0