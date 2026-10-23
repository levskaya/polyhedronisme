import { defineConfig } from 'vite';
import pegjs from 'pegjs';
import MarkdownIt from 'markdown-it';
import markdownItImsize from 'markdown-it-imsize';

const md = new MarkdownIt({
  html: true,
  breaks: true,
}).use(markdownItImsize);

function pegjsPlugin() {
  return {
    name: 'vite-plugin-pegjs',
    transform(source, id) {
      if (!id.endsWith('.pegjs')) return null;
      const compiled = pegjs.generate(source, {
        format: 'bare',
        output: 'source',
        optimize: 'speed',
      });
      return {
        code: `export default ${compiled};`,
        map: null,
      };
    },
  };
}

function markdownPlugin() {
  return {
    name: 'vite-plugin-markdown',
    transform(source, id) {
      if (!id.endsWith('.md')) return null;
      const html = md.render(source);
      return {
        code: `export default ${JSON.stringify(html)};`,
        map: null,
      };
    },
  };
}

export default defineConfig({
  publicDir: 'public',
  plugins: [pegjsPlugin(), markdownPlugin()],
  server: {
    port: 8080,
  },
});
