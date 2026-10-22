-- Keep each report image and its source caption in one non-floating box.
function Para(block)
  if #block.content == 1 and block.content[1].t == 'Image' then
    local caption = block.content[1].caption
    return {
      pandoc.RawBlock('latex', '\\par\\bigskip\\noindent\\begin{minipage}{\\linewidth}'),
      block,
      pandoc.Para(pandoc.Emph(caption)),
      pandoc.RawBlock('latex', '\\end{minipage}\\par\\bigskip')
    }
  end
end
