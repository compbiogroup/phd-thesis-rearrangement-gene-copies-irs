-- point links between doc pages to the PDFs
function Link(el)
  el.target = el.target:gsub("^([%w_]+)%.md$", "%1.pdf"):gsub("^%.%./BUGS%.md$", "../../BUGS.md")
  return el
end
-- "\|" escapes a pipe inside table cells on GitHub; pandoc keeps the backslash in code
function Code(el)
  el.text = el.text:gsub("\\|", "|")
  return el
end
