"""
Launch a live interactive Bonito.jl + WGLMakie web server for PlasmaWaves.jl.

Usage:
    julia --project=web web/serve.jl
"""

using Bonito
using WGLMakie
using PlasmaWaves

server = serve_website(; host = "127.0.0.1", port = 8080)
println("Press Ctrl+C to stop the PlasmaWaves.jl server.")
wait(server)
