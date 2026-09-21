--[[
  skip-gifs-in-pdf.lua

  LaTeX (used by Quarto for PDF rendering) does not support .gif images and
  will abort with "Unknown graphics extension: .gif". This filter runs only
  during PDF rendering and rewrites or drops .gif references so the PDF build
  succeeds.

  Behaviour:
    - For pattern animations (`.../animation.gif`), swap the src for the
      sibling `frame_0000.png`, which exists in every pattern directory.
    - For any other `.gif`, drop the image entirely (HTML still renders it;
      PDF just omits it).

  HTML rendering is unaffected — the filter is a no-op for non-PDF formats.
]]

function Image(img)
  if not quarto.doc.is_format("pdf") then
    return nil
  end

  local src = img.src
  if not src:match("%.gif$") then
    return nil
  end

  -- Pattern animations: swap to the static first-frame PNG sibling.
  if src:match("/animation%.gif$") then
    img.src = src:gsub("/animation%.gif$", "/frame_0000.png")
    return img
  end

  -- No known fallback — drop the image so LaTeX doesn't choke.
  return {}
end
