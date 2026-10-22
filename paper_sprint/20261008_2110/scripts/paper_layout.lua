-- Keep manuscript's explicit bilingual captions; prevent caption-only page tails.
function Blocks(blocks)
  local out = pandoc.List()
  local i = 1
  while i <= #blocks do
    local b = blocks[i]
    local s = pandoc.utils.stringify(b)
    if b.t == 'Header' and b.level == 2 and s:match('^5 ') then
      out:insert(pandoc.RawBlock('latex', '\\clearpage'))
    end
    if b.t == 'Para' and (s:match('^表[13]') or s:find('验证数值表', 1, true)) and blocks[i+1] and blocks[i+1].t == 'Table' then
      out:insert(pandoc.RawBlock('latex', '\\sprintneedspace{19\\baselineskip}'))
    end
    if b.t == 'Para' and #b.content == 1 and b.content[1].t == 'Image' and blocks[i+1] and blocks[i+1].t == 'Para' then
      out:insert(pandoc.RawBlock('latex', '\\begin{samepage}'))
      out:insert(b)
      out:insert(pandoc.RawBlock('latex', '\\nopagebreak[4]'))
      out:insert(blocks[i+1])
      out:insert(pandoc.RawBlock('latex', '\\end{samepage}'))
      i = i + 2
    else
      out:insert(b)
      i = i + 1
    end
  end
  return out
end
